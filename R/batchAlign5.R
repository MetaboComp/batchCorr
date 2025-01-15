#' BA: Aggregate presense/missingness per batch
#'
#' batchFlag aggregates presentness/missingness per feature in batches within 
#' the combination: batch x sample type
#'
#' @param PTnofill a peaktable with missing values 
#' (ie without either hard filling or imputation)
#' @param peakInfo a matrix with m/z (col1) and rt (col2) of features (rows)
#' @param NAhard proportion of NAs within batch for feature to be 
#' considered missing
#' @param batch vector (length = nSamples) containing batch information 
#' (e.g. A, B, C)
#' @param sampleGroup vector (length = nSamples) containing sample group 
#' information (e.g. QC, Sample, Reference)
#'
#' @return meta: metadata aggretated by batch x sample type
#' @return flagHard: a table of presentness/missingness per batch x sample type
#' @return mz: m/z values of the features
#' @return rt: rt values of the features
#' @noRd
batchFlag <- function(PTnofill,
                        batch,
                        sampleGroup,
                        peakInfo,
                        NAhard) {
    uniqBatch <- unique(batch)
    uniqGrp <- unique(sampleGroup)
    n <- length(uniqBatch) * length(uniqGrp)
    flagHard <- matrix(nrow = n, ncol = ncol(PTnofill))
    colnames(flagHard) <- colnames(PTnofill)
    batchMeta <- matrix(nrow = n, ncol = 2)
    colnames(batchMeta) <- c("batch", "sampleGroup")
    i <- 0
    for (b in uniqBatch) {
        for (g in uniqGrp) {
            i <- i + 1
            PTsub <- PTnofill[batch == b & sampleGroup == g,, drop = FALSE]
            NAs <- apply(PTsub, 2, function(x) sum(is.na(x)) / length(x))
            flagHard[i, ] <- ifelse(NAs >= NAhard, 1, 0)
            batchMeta[i, ] <- c(b, g)
        }
    }
    # Convert flags to !NA flags
    flagHard <- 1 - flagHard
    return(batchFlag <- list(
        meta = batchMeta,
        flagHard = flagHard,
        mz = peakInfo[, 1],
        rt = peakInfo[, 2]
    ))
}

#' BA: Find alignment candidates
#'
#' align will find candidates for alignment among features that have systematic 
#' missingness "within box" of maximum m/z & rt differences. Used internally.
#' @param flags a table of presentness/missingness per batch x sample type
#' @param mz m/z values of the features
#' @param rt rt values of the features
#' @param mzdiff maximum distance in m/z to be considered for alignment
#' @param rtdiff maximum distance in rt to be considered for alignment
#' @return events: all possible alignments
#' @return features: all features involved in possible alignments
#' @return clusters: groups of features linked through mutual events
#' @noRd
align <- function(flags, mz, rt, mzdiff, rtdiff) {
    # Take out misalignment candidates
    misIndex <- which(colSums(flags) != nrow(flags) & colSums(flags) != 0)
    event <- list()
    c <- 0
    for (i in seq_len(length(misIndex) - 1)) {
        # Take out next available misalignment candidate
        m <- misIndex[i]
        # Look at the ones with higher index
        misPlus <- misIndex[misIndex > m]
        # Their rts and mzs
        rtPlus <- rt[misPlus]
        mzPlus <- mz[misPlus]
        # Which are in the mz*rt box?
        inBox <- which(abs(rtPlus - rt[m]) < rtdiff &
                            abs(mzPlus - mz[m]) < mzdiff)

        if (length(inBox) > 0) {
            for (j in seq_along(inBox)) {
                n <- misPlus[inBox[j]]
                # Check overlapping features within box
                dotProd <- sum(flags[, m] * flags[, n])
                # For systematic misalignment to occur
                # these must be missed between features, ie dotProd=0
                if (dotProd == 0) {
                    c <- c + 1
                    event[[c]] <- c(
                        c, m, n, dotProd, sum(flags[, m]), sum(flags[, n]),
                        sum(flags[, m]) + sum(flags[, n])
                    ) # Take out data
                }
            }
        }
    }
    if (length(event) == 0) {
        return(list())
    }
    event <- matrix(unlist(event), ncol = 7, byrow = TRUE) # Rearrange
    event <- as.data.frame(event)
    colnames(event) <- c("event", "m", "n", "dotProd", "mSum", "nSum", "newSum")
    uniqFeat <- unique(c(event$m, event$n))
    allMN <- c(event$m, event$n)
    uniqAll <- unique(allMN)
    duplAll <- allMN[duplicated(c(allMN))]
    uniqM <- unique(event$m)
    duplM <- event$m[duplicated(event$m)]
    uniqN <- unique(event$n)
    duplN <- event$n[duplicated(event$n)]
    duplMN <- event$n[event$n %in% event$m]
    evTemp <- event
    features <- as.data.frame(matrix(nrow = 1, ncol = 8 + nrow(flags)))
    colnames(features) <- c(
        "cluster", "featureIndex", "mz", "rt", "to", "from", "MN",
        "NN", paste("bF", seq_len(nrow(flags)), sep = "")
    )
    ## Cluster together same m's
    for (j in seq_along(uniqM)) {
        uF <- uniqM[j]
        mNum <- sum(event$m == uF)
        nNum <- sum(event$n == uF)
        featIndex <- c(j, uF, mz[uF], rt[uF], 0, 0, 0, 0, flags[, uF])
        features <- rbind(features, featIndex)
        for (mIndex in seq_len(mNum)) {
            eventNum <- which(event$m == uF)[mIndex]
            nM <- event$n[eventNum]
            evTemp[eventNum, ] <- NA
            MN <- ifelse(nM %in% duplMN, 1, 0)
            NN <- ifelse(nM %in% duplN, 1, 0)
            featIndex <- c(j, nM, mz[nM], rt[nM], 0, 1, MN, NN, flags[, nM])
            features <- rbind(features, featIndex)
        }
    }
    features <- features[-1, ]
    ## cluster together same n's
    for (j in duplN) {
        clu <- unique(features[features$featureIndex %in% j, 1])
        for (k in 2:length(clu)) {
            features[features[, 1] %in% clu[k], 1] <- clu[1]
        }
    }
    # Add flag for combination of features
    features[features$featureIndex %in% duplN, 6] <-
        features[features$featureIndex %in% duplN, 6] + 1
    
    ## cluster together common m&n's
    for (j in duplMN) {
        clu <- unique(features[features$featureIndex %in% j, 1])
        for (k in 2:length(clu)) {
            features[features[, 1] %in% clu[k], 1] <- clu[1]
        }
    }
    # Add flag for combination of features
    features[features$featureIndex %in% duplMN, 5] <-
        features[features$featureIndex %in% duplMN, 5] + 1
    
    features <- features[!duplicated(features$featureIndex), ]
    features$cluster <- as.numeric(as.factor(features$cluster))
    clusters <- clusterMatrix(features)

    return(list(events = event, features = features, clusters = clusters))
}

#' BA: Sort features (alignment candidates) into clusters
#'
#' clusterMatrix will sort features into clusters. Used internally.
#' @param features all features involved in possible alignments
#' @return clusters: groups of features linked through mutual events
#' @noRd
clusterMatrix <- function(features) {
    nClust <- length(unique(features$cluster))
    clustMat <- 
        matrix(nrow = nClust, ncol = 3 + 
                length(grep("bF", colnames(features))))
    
    colnames(clustMat) <- c(
        "cluster", "nFeat", "dotProd",
        paste("bF", seq_along(grep("bF", colnames(features))), sep = "")
    )
    for (c in seq_len(nClust)) {
        clustSub <- features[features$cluster == c, ]
        vectMult <- as.numeric(clustSub[1, grep("bF", colnames(features))])
        for (f in 2:nrow(clustSub)) {
            dotProd <- sum(vectMult *
                            clustSub[f, grep("bF", colnames(features))])
            
            vectMult <- vectMult +
                as.numeric(clustSub[f, grep("bF", colnames(features))])
        }
        clustMat[c, ] <- c(c, nrow(clustSub), dotProd, vectMult)
    }
    return(as.data.frame(clustMat))
}

#' BA: Split overcrowded clusters into subclusters
#'
#' clustSplit will split clusters by iteratively removing maximum distances 
#' until "single possible candidates" are identified. The function works 
#' recursively of severely overcrowded clusters. Used internally.
#' @param clustFlags a table of presentness/missingness per batch x sample 
#' type for the features within cluster
#' @param mz m/z values of the features
#' @param rt rt values of the features
#' @return A vector with subcluster identifiers
#' @importFrom stats dist
#' @noRd
clustSplit <- function(clustFlags, mz, rt) {
    clustList <- rep(1, ncol(clustFlags))
    dotFlag <- matrix(0, ncol = ncol(clustFlags), nrow = ncol(clustFlags))
    for (i in seq_len(ncol(clustFlags))) {
        for (j in seq_len(ncol(clustFlags))) {
            dotFlag[i, j] <- ifelse(sum(clustFlags[, i] *
                                            clustFlags[, j]) == 0, 1, 0)
            
        }
    }
    rtDist <- as.matrix(dist(rt, diag = TRUE, upper = TRUE))
    mzDist <- as.matrix(dist(mz, diag = TRUE, upper = TRUE))
    heurFlag <- dotFlag * mzDist * rtDist
    heurFlag <- ifelse(heurFlag == 0, NA, heurFlag)
    hFlag <- heurFlag
    nDist <- apply(hFlag, 1, function(x) sum(!is.na(x)))
    while (min(nDist) > 1) {
        hFlag[which(hFlag == max(hFlag, na.rm = TRUE))] <- NA
        nDist <- apply(hFlag, 1, function(x) sum(!is.na(x)))
    }
    cand <- which(nDist == 1)
    coCand <- numeric()
    heur <- numeric()
    for (n in cand) {
        coCand <- c(coCand, which.min(hFlag[, n]))
        heur <- c(heur, min(hFlag[, n], na.rm = TRUE))
    }
    f1 <- cand[which.min(heur)]
    f2 <- coCand[which.min(heur)]
    clustList <- rep(2, ncol(clustFlags))
    clustList[cand[which.min(heur)]] <- clustList[coCand[which.min(heur)]] <- 1
    hFlagNew <- heurFlag[
        !seq_len(nrow(heurFlag)) %in% c(f1, f2),
        !seq_len(nrow(heurFlag)) %in% c(f1, f2)
    ]
    if (is.matrix(hFlagNew)) {
        nDistNew <- apply(hFlagNew, 1, function(x) sum(!is.na(x)))
        if (sum(nDistNew) == 0) {
            clustList[clustList == 2] <- seq(seq_len(sum(clustList == 2))) + 1
        } else {
            subFeats <- as.data.frame(cbind(clustList, t(clustFlags)))
            colnames(subFeats)[1] <- "cluster"
            subClustMat <- clusterMatrix(subFeats)
            if (subClustMat$dotProd[2] != 0) {
                ### clustSplit(new matrix)
                message("Second round of clustSplit highly experimental")
                newSplit <- clustSplit(
                    clustFlags[, clustList == 2],
                    mz[clustList == 2],
                    rt[clustList == 2]
                )
                clustList[clustList == 2] <- newSplit + 1
            }
        }
    }
    return(clustList)
}

#' BA: Plot the identified clusters
#'
#' plotClust will plot clusters identified based on sample type. Used 
#' internally.
#' @param batchflag a table of presentness/missingness per batch x sample 
#' type for the features within cluster
#' @param grpFlag a vector corresponding to the rows of batchFlag used for 
#' sample type
#' @param cluster extracted cluster to be plotted
#' @param text cluster identifier
#' @param color vector with colors
#' @param mzwidth plot span of m/z
#' @param rtwidth plot span of rt
#' @importFrom graphics points
#' @noRd
plotClust <- function(batchflag,
                        grpFlag,
                        cluster,
                        text,
                        color = 2,
                        mzwidth = 0.02,
                        rtwidth = 100) {
    bF <- batchflag
    mzspan <- c(mean(cluster$mz) - mzwidth / 2, mean(cluster$mz) + mzwidth / 2)
    rtspan <- c(mean(cluster$rt) - rtwidth / 2, mean(cluster$rt) + rtwidth / 2)
    plot(cluster$rt,
        cluster$mz,
        col = color,
        pch = 16,
        xlim = rtspan,
        ylim = mzspan,
        xlab = "Retention time (s)",
        ylab = "m/z"
    )
    text(rtspan[1],
        mzspan[2] - mzwidth / 20,
        text,
        cex = 1.3,
        pos = 4
    )
    bFs <- apply(
        cluster[, grep("bF", colnames(cluster))], 1,
        function(x) paste(x, collapse = " ")
    )
    text(cluster$rt,
        cluster$mz,
        cluster$featureIndex,
        pos = 3
    )
    text(cluster$rt,
        cluster$mz,
        bFs,
        pos = 1
    )

    whichPlot <- (abs(bF$mz - mean(cluster$mz)) < mzwidth / 2 &
        abs(bF$rt - mean(cluster$rt)) < rtwidth / 2)
    mzPlot <- bF$mz[whichPlot]
    rtPlot <- bF$rt[whichPlot]
    bfPlot <- apply(
        bF$flagHard[grpFlag, whichPlot], 2,
        function(x) paste(x, collapse = " ")
    )

    points(rtPlot,
        mzPlot,
        col = 1,
        pch = 1
    )
    text(rtPlot,
        mzPlot,
        which(whichPlot),
        pos = 3
    )
    text(rtPlot,
        mzPlot,
        bfPlot,
        pos = 1
    )
}

#' BA: Find alignment candidates
#'
#' alignIndex will find features systematically misaligned between batches 
#' using "sample type" information.
#' @param batchflag a table of presentness/missingness per batch x 
#' sample type for the features within cluster
#' @param grpType sample type to be used to find batch alignments
#' @param mzdiff maximum distance in m/z to be considered for alignment
#' @param rtdiff maximum distance in rt to be considered for alignment
#' @param report boolean for plotting results of owercrowded cluster 
#' splitting for visual examination
#' @param reportName string with name of report file
#' @return An object (list) consisting of the following features:
#' @return grpType: sample type used to find batch alignments
#' @return shift: list with features to combine
#' @return events: all possible alignment candidates (events)
#' @return features: features after splitting overcrowded clusters
#' @return clusters: clusters after splitting overcrowded clusters
#' @return oldFeatures: features before splitting overcrowded clusters
#' @return oldClusters: clusters before splitting overcrowded clusters
#' @importFrom grDevices dev.off pdf
#' @noRd
alignIndex <- function(batchflag,
                        grpType,
                        mzdiff,
                        rtdiff,
                        report,
                        reportName = "cluster_splits",
                        reportPath) {
    bF <- batchflag
    if (!grpType %in% bF$meta[, "sampleGroup"]) {
        stop("Group ", grpType, " not in metadata.")
    }
    # Take out matrix based on group type (such as QC or long-term Ref)
    grpSub <- bF$flagHard[bF$meta[, 2] == grpType, ]
    mz <- bF$mz
    rt <- bF$rt
    a2 <- a1 <- align(grpSub, mz, rt, mzdiff = mzdiff, rtdiff = rtdiff)
    a1Clust <- a1$clusters
    if (is.null(a1Clust)) {
        stop("There are no alignment candidates.
        Therefore, between-batch alignment is not possible.
        Consider expanding mzdiff and/or rtdiff.")
    }
    splits <- which(a1Clust$dotProd != 0)
    if (report) pdf(file = paste(reportPath, reportName, ".pdf", sep = ""))
    for (s in splits) {
        cluster <- a1$features[a1$features$cluster == s, ]
        # print(cluster)
        mzSp <- cluster$mz
        rtSp <- cluster$rt
        flagSp <- t(cluster[9:ncol(a1$features)])
        alignSplit <- clustSplit(flagSp, mzSp, rtSp)
        text <- paste("Original cluster", s)
        if (report) {
            plotClust(bF,
                        bF$meta[, 2] == grpType,
                        cluster,
                        text = text,
                        color = alignSplit + 1)
        }
        newClust <- ifelse(alignSplit == 1,
                            s,
                            max(a2$features$cluster) + alignSplit - 1)
        
        a2$features$cluster[a2$features$cluster == s] <- newClust
    }
    if (report) dev.off()
    ### Remove duplicates
    a2$features <-
        a2$features[a2$features$cluster %in% 
                unique(a2$features$cluster[duplicated(a2$features$cluster)]), ]
    
    a2$features <- a2$features[order(a2$features$cluster), ]
    a2$features$cluster <- as.numeric(as.factor(a2$features$cluster))
    a2$clusters <- clusterMatrix(a2$features)
    ### Make index string for alignment
    nClust <- nrow(a2$clusters)
    shiftList <- numeric(ncol(bF$flagHard))
    shiftGrp <- character(ncol(bF$flagHard))
    for (c in seq_len(nClust)) {
        feats <- a2$features$featureIndex[a2$features$cluster == c]
        shiftList[feats] <- feats[1]
        shiftGrp[feats] <- as.character(grpType)
        shift <- data.frame(list = shiftList, sampleGroup = shiftGrp)
    }
    return(list(
        grpType = grpType,
        shift = shift,
        events = a2$events,
        features = a2$features,
        clusters = a2$clusters,
        oldFeatures = a1$features,
        oldClusters = a1$clusters
    ))
}

#' BA: Plot clusters of aligned features
#'
#' plotAlign will plot clusters of aligned features
#' @param batchflag a table of presentness/missingness per batch x sample 
#' type for the features within cluster
#' @param alignindex An object (list) consisting of alignment information
#' @param clust which cluster(s) to plot. If missing, plots all clusters.
#' @param plotType whether to plot to list ("plot") or to pdf ("pdf") 
#' Default: "plot"
#' @param reportName name of plotfile
#' @param mzwidth plot span of m/z
#' @param rtwidth plot span of rt
#' @noRd
plotAlign <- function(batchflag,
                        alignindex,
                        clust,
                        reportName = "aligned_clusters",
                        reportPath) {
    aI <- alignindex
    bF <- batchflag
    pdf(file = paste(reportPath, reportName, ".pdf", sep = ""))
    if (missing(clust)) {
        clustPlots <- seq_len(dim(aI$clusters)[1])
    } else {
        clustPlots <- clust
    }
    # cat(clustPlots)
    if (length(clustPlots) == 2) {
        par(mfrow = c(2, 1))
    }
    if (length(clustPlots) > 2) {
        par(mfrow = c(3, 2))
    }
    for (c in clustPlots) {
        cluster <- aI$features[aI$features$cluster == c, ]
        # print(cluster)
        mz <- cluster$mz
        rt <- cluster$rt
        flags <- t(cluster[9:ncol(aI$features)])
        text <- paste("Cluster", c)
        plotClust(
            batchflag = bF,
            grpFlag = bF$meta[, 2] == aI$grpType,
            cluster = cluster,
            text = text,
            color = 2
        )
    }
    dev.off()
}

#' BA: Alignment of peaktable based on alignIndex and batchFlag data
#'
#' batchAlign will use the 'batchFlag' and 'alignIndex' information to align 
#' sample peaks that are systematically misaligned across batches.
#'
#' @param batchflag a table of presentness/missingness per batch x sample type 
#' for the features within cluster
#' @param alignindex An object (list) consisting of alignment information
#' @param peaktable_filled a peaktable without missing values 
#' (ie after hard filling or imputation)
#' @param batch vector (length = nSamples from peaktable_filled) 
#' containing batch information (e.g. A, B, C)
#'
#' @return An object (list) consisting of the following features:
#' @return PTalign: A peaktable with batch aligned peak areas
#' @return boolAveragedAlign: boolean vector of features where alignment has 
#' been made using feature averaging 
#' (i.e. where batches are missing within features). 
#' Length: same as final number of features
#' @return PTfill: peaktable without missing values (indata)
#' @return boolKeep: boolean vector of features kept after alignment 
#' (the rest of the combined features are deleted)
#' @return boolAveragedFill: boolean vector of features where alignment has been
#'  made using feature averaging 
#'  (i.e. where batches are missing within features). 
#'  Length: same as original number of features
#' @return aI: alignIndex object (indata)
#' @importFrom stats aggregate
#' @noRd
batchAlign <- function(batchflag,
                        alignindex,
                        peaktable_filled,
                        batch) {
    bF <- batchflag
    flags <- bF$flagHard
    batchFlags <- as.matrix(aggregate(flags, list(bF$meta[, 1]), "sum")[, -1])
    bFbatch <- bF$meta[, 1]
    uniqBatch <- unique(bFbatch)
    aI <- alignindex
    shift <- aI$shift$list
    uniqSh <- unique(shift[shift != 0])
    grp <- aI$shift$grp
    PTfill <- peaktable_filled
    boolKeep <- rep(TRUE, ncol(PTfill))
    boolAveraged <- !boolKeep
    for (sh in uniqSh) {
        newFeat <- numeric(nrow(PTfill))
        feats <- which(shift == sh)
        subVect <- flags[, feats]
        if (is.matrix(subVect)) {
            vectAdd <- subVect[, 1]
            subBatchVect <- batchFlags[, feats]
            batchVectAdd <- subBatchVect[, 1]
            subFeats <- PTfill[, feats]
            aveFeat <- 
                rowSums(t(colSums(subVect) * t(subFeats))) / sum(subVect)
            for (f2 in 2:length(feats)) {
                vect2 <- subVect[, f2]
                dotProd <- sum(vectAdd * vect2, na.rm = TRUE)
                vectAdd <- vectAdd + vect2
                batchVect2 <- subBatchVect[, f2]
                batchDotProd <- sum(batchVectAdd * batchVect2, na.rm = TRUE)
                batchVectAdd <- batchVectAdd + batchVect2
            }
            if (dotProd == 0 & batchDotProd == 0) {
                for (f in feats) {
                    boolFeat <- rep(FALSE, nrow(PTfill))
                    bFlag <- which(batchFlags[, f] > 0)
                    for (bFl in bFlag) {
                        batchID <- uniqBatch[bFl]
                        boolFeat <- 
                            boolFeat | ifelse(batch == batchID, TRUE, FALSE)
                    }
                    newFeat[boolFeat] <- PTfill[boolFeat, f]
                }
                ### How to deal with non-present batches???
                zeros <- which(rowSums(subBatchVect) == 0)
                if (length(zeros) > 0) {
                    boolAveraged[feats] <- TRUE
                    boolFeat <- rep(FALSE, nrow(PTfill))
                    for (z in zeros) {
                        batchID <- uniqBatch[z]
                        boolFeat <- 
                            boolFeat | ifelse(batch == batchID, TRUE, FALSE)
                    }
                    newFeat[boolFeat] <- aveFeat[boolFeat]
                }
                PTfill[, feats[1]] <- newFeat
                PTfill[, feats[-1]] <- 0
                boolKeep[feats[-1]] <- FALSE
            }
        }
    }
    PTfill <- PTfill[, boolKeep]
    return(list(
        PTalign = PTfill,
        boolAveragedAlign = boolAveraged[boolKeep],
        PTfill = peaktable_filled,
        boolKeep = boolKeep,
        boolAveragedFill = boolAveraged,
        aI = alignindex
    ))
}
