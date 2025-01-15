#' Coefficient of variation (CV)
#'
#' Calculates CV per column in a matrix. Aka relative standard deviation (RSD).
#' @param mat a matrix with variables as columns
#' @return a vector of CVs
#' @noRd
## Simple function for calculating cv per column (ie variable)
cv <- function(mat) {
    if (is.null(dim(mat))) {
        cv <- abs(sd(mat) / mean(mat))
    } else {
        mean <- apply(mat, 2, mean)
        sd <- apply(mat, 2, sd)
        cv <- abs(sd / mean)
    }
    return(cv)
}

#' Root mean squared distance
#'
#' Calculate root mean squared distance from center point
#' @param mat a matrix containing observations as rows and variables as columns
#' @return a numeric rmsDist
#' @noRd
## Simple function for calculating root mean square distance
rmsDist <- function(mat) {
    mean <- colMeans(mat)
    rmsd <- sqrt(sum(apply(mat, 1, function(x) sum((x - mean)^2))) / nrow(mat))
    return(rmsd)
}

#' DC: Cluster features
#'
#' clust will perform clustering of features with similar drift pattern by 
#' projecting scaled QC features as coordinates in observation/injection space
#' @param QCInjs vectors with QC injections in sequence
#' @param QCFeats Feature matrix for QC injections: Scaled (but not centered) 
#' features as columns; Injections as rows
#' @param modelNames Which 'mclust' models to consider 
#' (see mclust package for details)
#' @param G Which number of clusters to consider 
#' (see mclust package for details)
#' @param report boolean whether to print a pdf report of clustering results
#' @return a clust object containing:
#' @return QCInjs: as indata
#' @return QCFeats: as indata
#' @return BIC: Bayesian Information Criteria-object for the 
#' investigated clustering models
#' @return clust: summary of the BIC-object (i.e. clustering)
#' @return BICTime: Time required for the BIC-calculation
#' @return clustTime: Time required for the clustering
#' @import mclust
#' @noRd
clust <- function(QCInjs,
                    QCFeats,
                    modelNames,
                    G,
                    report,
                    reportPath) {
    if (length(QCInjs) != nrow(QCFeats)) {
        stop("nrow(QCFeats) not equal to length(QCInjs)")
    }
    message("Mclust ")
    combinations <- expand.grid("modelNames" = modelNames, "G" = G)
    startTime <- proc.time()[3]
    merged_models <- NULL
    models <- BiocParallel::bplapply(
        seq_len(nrow(combinations)),
        function(comb, combinations, QCFeats, merged_models) {
            mclBIC <- mclustBIC(t(QCFeats),
                G = combinations$G[comb],
                modelNames = combinations$modelNames[comb]
            )
        },
        combinations, QCFeats, merged_models
    )
    endTime <- proc.time()[3]
    BICtime <- endTime - startTime
    merged_models <- NULL
    # In case of missing BIC values for models, the failed merging is handled
    for (i in seq_along(models)) {
        tryCatch(
            {
                merged_models <- mclustBICupdate(merged_models, models[[i]])
            },
            error = function(e) {
                if (grepl("missing value where TRUE/FALSE needed", e)) {
                    message(
                "Missing BIC values for some model(s). You may want to check",
                " which models have missing BIC values and/or try again."
                    )
                }
            }
        )
    }
    startTime <- proc.time()[3]
    MC <- summary(merged_models, data = t(QCFeats))
    endTime <- proc.time()[3]
    clustTime <- endTime - startTime
    if (report == TRUE) {
        pdf(file = paste(reportPath,
            "cluster_BIC_",
            format(Sys.time(), format = "%y%m%d_%H%M"),
            ".pdf",
            sep = ""
        ))
        plot(merged_models)
        dev.off()
    }
    message(
        "MClust final model with ", MC$G,
        " clusters and ", MC$modelName, " geometry."
    )
    message(
        "BIC performed in ", BICtime,
        " seconds and clustering in ", clustTime, " seconds."
    )
    return(list(
        QCInjs = QCInjs,
        QCFeats = QCFeats,
        BIC = merged_models,
        clust = MC,
        BICTime = BICtime,
        clustTime = clustTime
    ))
}
.calc_driftCalc <- function(nclass,
                            classes,
                            varClust,
                            QCFeats,
                            QCInjs,
                            injs,
                            spar,
                            corMat,
                            deltaDist,
                            rmsdRaw,
                            cvRaw,
                            cvs,
                            ratios,
                            cvCorr,
                            smoothFunc,
                            report,
                            reportPath) {
    rmsdRaw <- rmsDist(QCFeats)
    for (n in seq_len(nclass)) {
        # Allocate matrix (for drift corrected variables) for later QC dist calc
        QCFeatsCorr <- QCFeats
        whichVars <- which(classes == n)
        vars <- QCFeats[, whichVars] # Take out cluster variables
        varClust[[n]] <- colnames(QCFeats)[whichVars]
        V <- as.data.frame(cbind(QCInjs, vars)) # Arrange and rearrange data
        V <- melt(V, id.vars = "QCInjs")
        V <- V[order(V$QCInjs), ]
        # Interpolations
        if (length(QCInjs) <= 3) { # 2nd degree polynomial if <4 data points
            Fit <- lm(value ~ poly(QCInjs, 2), data = V)
            Pred <- predict(Fit, data.frame(QCInjs = injs))
            Pred <- data.frame(x = injs, y = Pred)
        } else {
            if (smoothFunc == "spline") {
                # Cubic spline regression otherwise
                splineFit <- smooth.spline(V$QCInjs, V$value, spar = spar)
                # Predict drift over all injections
                Pred <- predict(splineFit, injs)
            } else {
                loessFit <- loess(value ~ QCInjs, data = V, span = spar)
                Pred <- predict(loessFit, data.frame(QCInjs = injs))
                Pred <- data.frame(x = injs, y = Pred)
            }
        }
        # Calculate correction factors for all injections
        corFact <- Pred$y[1] / Pred$y
        # Store cluster correction factors in "master" matrix
        corMat[, n] <- corFact
        # Bring out correction factors for QC samples specifically
        corQC <- corFact[QCInjs - min(QCInjs) + 1]
        # correct drift within cluster
        QCFeatsCorr[, classes == n] <- QCFeats[, classes == n] * corQC
        ## Calculate rmsDist
        rmsdCorr <- rmsDist(QCFeatsCorr)
        deltaDist[n] <- rmsdCorr - rmsdRaw

        # Calculate change in average QC distance
        # deltaDist[n]=mean(dist(QCFeatsCorr))-meanDistQCFeats

        cvRaw[n] <- mean(cv(QCFeats[, classes == n]))
        cvCorr[n] <- mean(cv(QCFeatsCorr[, classes == n]))
        cvs[[n]] <- data.frame(
            Raw = cv(QCFeats[, classes == n]),
            Corr = cv(QCFeatsCorr[, classes == n])
        )
        ratios[n, ] <- c(
            sum(cvs[[n]]$Raw < 0.15) / nrow(cvs[[n]]),
            sum(cvs[[n]]$Corr < 0.15) / nrow(cvs[[n]]),
            sum(cvs[[n]]$Raw < 0.2) / nrow(cvs[[n]]),
            sum(cvs[[n]]$Corr < 0.2) / nrow(cvs[[n]])
        )
        if (report == TRUE) {
            # Plot drift and drift function
            matplot(QCInjs,
                QCFeats[, classes == n],
                type = "l",
                lty = 1,
                col = "grey",
                ylim = range(QCFeats[, classes == n]),
                main = paste("Cluster ", n, "; n=", sum(classes == n),
                    "; Raw; Mean CV=",
                    round(mean(cv(QCFeats[, classes == n])), 3),
                    sep = ""
                ),
                ylab = "Scaled intensity",
                xlab = "Injection number"
            )
            lines(Pred, pch = 2)
            matplot(QCInjs, QCFeatsCorr[, classes == n],
                type = "l",
                lty = 1, col = "grey", ylim = range(QCFeats[, classes == n]),
                main = paste("Corrected; Mean CV=",
                    round(mean(cv(QCFeatsCorr[, classes == n])), 3),
                    sep = ""
                ),
                ylab = "Scaled intensity",
                xlab = "Injection number"
            )
        }
    }
    return(list(
        "ratios" = ratios,
        "cvRaw" = cvRaw,
        "cvCorr" = cvCorr,
        "corMat" = corMat,
        "deltaDist" = deltaDist,
        "varClust" = varClust
    ))
}

#' DC: Calculate intensity drift per cluster
#'
#' Clustered, scaled QC features are pooled together and intensity 
#' drift patterns are modelled per cluster
#' @param QCClust a clust object
#' @param smoothFunc choice of regression function: spline or loess 
#' (defaults to spline)
#' @param spar smoothing parameter for spline or loess (defaults to 0.2)
#' @param report boolean whether to print a pdf report of drift models
#' @return a driftCalc object containing:
#' @return original information from clust object (indata) and
#' @return actionInfo
#' @return ratios
#' @return corMat
#' @return deltaDist
#' @return varClust
#' @importFrom stats lm loess smooth.spline predict
#' @importFrom graphics lines matplot par
#' @import reshape
#' @noRd
driftCalc <- function(QCClust,
                        smoothFunc = c("spline", "loess"),
                        spar,
                        report,
                        reportPath) {
    smoothFunc <- match.arg(smoothFunc)
    MC <- QCClust$clust
    QCInjs <- QCClust$QCInjs
    QCFeats <- QCClust$QCFeats
    if (length(QCInjs) != nrow(QCFeats)) {
        stop("nrow(QCFeats) not equal to length(QCInjs)")
    }
    # Extract classes
    
    # Take out total number of identified
    # clusters/components/groups/classes/whatever you want to call them
    nclass <- MC$G 
    
    # Take out the classifications for the different variables
    classes <- MC$classification
    
    # Allocate variables
    # allocate vector for effect size of drift correction per cluster
    cvRaw <- cvCorr <- deltaDist <- numeric(nclass) 
    # Make injection list
    injs <- min(QCInjs):max(QCInjs) 
    # Allocate matrix with correction function (rows) per cluster (column)
    corMat <- matrix(nrow = length(injs), ncol = nclass) 
    cvs <- varClust <- list()
    ratios <- matrix(nrow = nclass, ncol = 4)
    rownames(ratios) <- paste("cluster", seq_len(nclass), sep = "")
    colnames(ratios) <- c("raw.15", "corr.15", "raw.2", "corr.2")
    if (report == TRUE) {
        pdf(file = paste(reportPath,
            "cluster_G",
            nclass, "_",
            format(Sys.time(), format = "%y%m%d_%H%M"),
            ".pdf",
            sep = ""
        ))
        par(mfrow = c(2, 1))
        par(mar = c(2, 4, 2, 0))
    }
    ### Calculate drift correction for each cluster
    # Calculate distance on scaled variables
    rmsdRaw <- rmsDist(QCFeats)
    driftCalc_list <- .calc_driftCalc(
        nclass = nclass, classes = classes,
        varClust = varClust, QCFeats = QCFeats,
        QCInjs = QCInjs, injs = injs, spar = spar,
        corMat = corMat, deltaDist = deltaDist,
        rmsdRaw = rmsdRaw, cvRaw = cvRaw,
        cvs = cvs, ratios = ratios, cvCorr = cvCorr,
        smoothFunc, report = report,
        reportPath = reportPath
    )

    if (report == TRUE) {
        dev.off() # Close pdf file
    }
    clustComm <- rep("None", nclass)
    actionInfo <- data.frame(
        number = seq_len(nclass),
        n = vapply(driftCalc_list$varClust, length, integer(1)),
        action = clustComm,
        CVRaw = driftCalc_list$cvRaw,
        CVCorr = driftCalc_list$cvCorr
    )
    QCClust$actionInfo <- actionInfo
    QCClust$ratios <- driftCalc_list$ratios
    QCClust$corMat <- driftCalc_list$corMat
    QCClust$deltaDist <- driftCalc_list$deltaDist
    QCClust$varClust <- driftCalc_list$varClust
    QCDriftCalc <- QCClust
    message("Calculation of QC drift profiles performed.")
    return(QCDriftCalc)
}

#' DC: Correct for intensity drift per cluster
#'
#' Perform signal intensity drift correction if resulting in increased 
#' quality of data 
#' (measured by reduced root-mean-squared distances of reference samples).
#' @param QCDriftCalc a DriftCalc object
#' @param refList a reference sample object with 'inj' and 'Feats'
#' @param refType at present, the options "one" and "none" are supported 
#' for one or no reference samples present.
#' @param CorrObj a batch object to be corrected with 'inj' and 'Feats'. 
#' If not present defaults to the QC object
#' @param report boolean whether to print a pdf report of drift models
#' @return A Corr object consisting of:
#' @return original information from DriftCalc object (indata) and
#' @return updated actionInfo (actions taken per cluster)
#' @return QCFeatsCorr Drift-corrected QC features
#' @return refType (indata)
#' @return drift corrected reference samples
#' @return drift corrected corr/batch samples
#' @importFrom graphics hist legend
#' @importFrom grDevices rgb
#' @noRd
## Perform drift correction for clusters IF rmsdRef is improved
driftCorr <- function(QCDriftCalc,
                        refList = NULL,
                        refType = c("none", "one"),
                        CorrObj = NULL,
                        report,
                        reportPath) {
    if (missing(refType)) refType <- "none"
    deltaDist <- QCDriftCalc$deltaDist
    varClust <- QCDriftCalc$varClust
    removeFeats <- QCDriftCalc$removeFeats
    if (!is.null(removeFeats)) {
        keepClust <- QCDriftCalc$keepClust
        corrQCTemp <- corrQC <- QCDriftCalc$QCFeatsClean
    } else {
        keepClust <- seq_len(QCDriftCalc$clust$G)
        corrQCTemp <- corrQC <- QCDriftCalc$QCFeats
    }
    injQC <- QCDriftCalc$QCInjs
    corMat <- QCDriftCalc$corMat
    clustComm <- as.character(QCDriftCalc$actionInfo$action)
    ordDist <- order(deltaDist)
    ordDist <- ordDist[ordDist %in% keepClust]
    if (refType == "one") {
        refClean <- refList$Feats[, !colnames(refList$Feats) %in% removeFeats]
        injRef <- refList$inj
        corrRefTemp <- corrRef <- refClean
        if (length(injRef) != nrow(corrRef)) {
            stop("nrow(corrRef) not equal to length(injRef)")
        }
    }
    if (missing(CorrObj)) {
        injTest <- injQC
        corrTest <- corrQC
        CorrObj <- list(inj = injTest, Feats = corrTest)
    } else {
        injTest <- CorrObj$inj
        corrTest <- CorrObj$Feats
        corrTest <- corrTest[, !colnames(corrTest) %in% removeFeats]
    }
    if (length(injTest) != nrow(corrTest)) {
        stop("nrow(corrTest) not equal to length(injTest)")
    }
    for (i in seq_along(keepClust)) {
        n <- ordDist[i]
        # take out cluster correction factors from "master" matrix
        corFact <- corMat[, n] 
        corrFeats <- varClust[[n]]
        # Scheme for (suboptimal) situation without Ref samples
        if (refType == "none") { 
            # Bring out correction factors for QC samples specifically
            corQC <- corFact[injQC - min(injQC) + 1] 
            corrQCTemp[, colnames(corrQC) %in% corrFeats] <- 
                corrQC[, colnames(corrQC) %in% corrFeats] * corQC
            
            if (rmsDist(corrQCTemp) < rmsDist(corrQC)) {
                clustComm[n] <- "Corr_QC"
                corrQC <- corrQCTemp
                # Bring out correction factors for QC samples specifically
                corQC <- corFact[injQC - min(injQC) + 1] 
                
                corrQC[, colnames(corrQC) %in% corrFeats] <- 
                    corrQC[, colnames(corrQC) %in% corrFeats] * corQC
                # Bring out correction factors for Test samples specifically
                corTest <- corFact[injTest - min(injQC) + 1] 
                corrTest[, colnames(corrTest) %in% corrFeats] <- 
                    corrTest[, colnames(corrTest) %in% corrFeats] * corTest
                
            } else {
                corrQCTemp <- corrQC
            }
        }
        # For situation with multiple injections of singular non-QC Ref samples
        if (refType == "one") { 
            corRef <- corFact[injRef - min(injQC) + 1]
            corrRefTemp[, colnames(corrRefTemp) %in% corrFeats] <- 
                corrRefTemp[, colnames(corrRefTemp) %in% corrFeats] * corRef
            
            if (rmsDist(corrRefTemp) < rmsDist(corrRef)) {
                clustComm[n] <- "Corr_1Ref"
                corrRef <- corrRefTemp
                # Bring out correction factors for QC samples specifically
                corQC <- corFact[injQC - min(injQC) + 1] 
                corrQC[, colnames(corrQC) %in% corrFeats] <- 
                    corrQC[, colnames(corrQC) %in% corrFeats] * corQC
                
                # Bring out correction factors for Test samples specifically
                corTest <- corFact[injTest - min(injQC) + 1] 
                corrTest[, colnames(corrTest) %in% corrFeats] <- 
                    corrTest[, colnames(corrTest) %in% corrFeats] * corTest
                
            } else {
                corrRefTemp <- corrRef
            }
        }
    }
    if (report == TRUE) {
        pdf(file <- paste(reportPath,
            "Hist_Corrected_",
            format(Sys.time(), format = "%y%m%d_%H%M"),
            ".pdf",
            sep = ""
        ))
        if (!is.null(removeFeats)) {
            cvBefore <- cv(QCDriftCalc$QCFeats)
        } else {
            cvBefore <- cv(QCDriftCalc$QCFeats)
        }
        histCombined <- hist(c(cvBefore, cv(corrQC)), 30, plot = FALSE)
        breaks <- histCombined$breaks
        histBefore <- hist(cvBefore, plot = FALSE, breaks = breaks)
        histAfter <- hist(cv(corrQC), plot = FALSE, breaks = breaks)
        ymax <- max(c(histBefore$counts, histAfter$counts))
        hist(cvBefore,
            breaks = breaks,
            ylim = c(0, ymax),
            col = rgb(0, 0, 0, 1),
            main = "Cluster correction",
            xlab = "CV (feature)"
        )
        hist(cv(corrQC),
            ylim = c(0, ymax),
            breaks = breaks,
            col = rgb(1, 1, 1, .5),
            add = TRUE
        )
        legend("topright",
            legend = c("Clean", "Corrected"),
            fill = c(rgb(0, 0, 0, 1), rgb(1, 1, 1, 0.5)), bty = "n"
        )
        dev.off()
    }
    QCDriftCalc$actionInfo$action <- clustComm
    QCDriftCalc$QCFeatsCorr <- corrQC
    QCDriftCalc$RefType <- refType
    if (refType == "one") {
        QCDriftCalc$RefInjs <- injRef
        QCDriftCalc$RefFeats <- refList$Feats
        QCDriftCalc$RefFeatsClean <- refClean
        QCDriftCalc$RefFeatsCorr <- corrRef
    }
    QCDriftCalc$TestInjs <- injTest
    QCDriftCalc$TestFeats <- CorrObj$Feats
    # QCDriftCalc$TestFeatsClean=
    QCDriftCalc$TestFeatsCorr <- corrTest
    QCCorr <- QCDriftCalc
    message(
        "\nDrift correction of ",
        sum(QCCorr$actionInfo$action != "None"), " out of ",
        QCCorr$clust$G, " clusters"
    )
    if (refType == "none") {
        message("using QC samples only.")
    } else {
        ("validated by external reference samples.")
    }
    message("\nCorrected peak table in $TestFeatsCorr\n")
    return(QCCorr)
}

#' DC: Remove features not passing QC test
#'
#' Remove those features with CV(QC)>limit (defined by user)
#' @param QCCorr a Corr object
#' @param CVlimit user-defined QC riterion to pass. Defaults to 0.2
#' @param report boolean whether to print a pdf report of drift models
#' @return A Clean object consisting of:
#' @return original information from Corr object (indata) and
#' @return updated actionInfo (actions taken per cluster)
#' @return QCFeatsFinal: Final peaktable (scaled) for QC samples
#' @return TestFeatsFinal: Final peaktable for batch samples
#' @return finalVars: All features kept within the final peaktable
#' @return QCcvs: CVs per features kept within the final peaktable
#' @noRd
## Remove individual variables with CV>limit
cleanVar <- function(QCCorr,
                        CVlimit,
                        report = FALSE,
                        reportPath) {
    QCFeats <- QCCorr$QCFeats
    removeFeats <- QCCorr$removeFeats
    if (!is.null(removeFeats)) {
        QCFeatsClean <- QCCorr$QCFeatsClean
    } else {
        QCFeatsClean <- QCCorr$QCFeats
    }
    QCFeatsCorr <- QCCorr$QCFeatsCorr
    cvIndex <- which(cv(QCFeatsCorr) > CVlimit)
    finalIndex <- rep(TRUE, ncol(QCFeatsCorr))
    if (length(cvIndex) > 0) finalIndex[cvIndex] <- FALSE
    QCFeatsFinal <- QCFeatsCorr[, finalIndex]
    if (QCCorr$RefType == "one") {
        RefFeatsFinal <- QCCorr$RefFeatsCorr[, finalIndex]
    }
    TestFeatsFinal <- QCCorr$TestFeatsCorr[, finalIndex]
    finalVars <- colnames(QCFeatsFinal)
    cvFeats <- mean(cv(QCFeats)) # 0.2276
    cvFeatsClean <- mean(cv(QCFeatsClean)) # 0.1455
    cvFeatsCorr <- mean(cv(QCFeatsCorr)) # 0.1144
    cvFeatsFinal <- mean(cv(QCFeatsFinal)) # 0.1070
    QCcvs <- data.frame(
        cvFeats = cvFeats,
        cvFeatsClean = cvFeatsClean,
        cvFeatsCorr = cvFeatsCorr,
        cvFeatsFinal = cvFeatsFinal
    )
    if (report == TRUE) {
        pdf(file = paste(reportPath,
            "Hist_Final_",
            format(Sys.time(), format = "%y%m%d_%H%M"),
            ".pdf",
            sep = ""
        ))
        histCombined <- hist(c(cv(QCFeatsClean), cv(QCFeatsFinal)),
                                30,
                                plot = FALSE)
        breaks <- histCombined$breaks
        histBefore <- hist(cv(QCFeatsClean), breaks = breaks, plot = FALSE)
        histAfter <- hist(cv(QCFeatsFinal), breaks = breaks, plot = FALSE)
        ymax <- max(c(histBefore$counts, histAfter$counts))
        hist(cv(QCFeatsClean),
            breaks = breaks,
            ylim = c(0, ymax),
            col = rgb(0, 0, 0, .75),
            main = "Cluster cleanup",
            xlab = "CV (feature)"
        )
        hist(cv(QCFeatsFinal),
            breaks = breaks,
            ylim = c(0, ymax),
            col = rgb(1, 1, 1, .8),
            add = TRUE
        )
        legend("topright",
            legend = c("Clean", "Final"),
            fill = c(rgb(0, 0, 0, 1), rgb(1, 1, 1, 0.5)),
            bty = "n"
        )
        dev.off()
    }
    if (QCCorr$RefType == "one") {
        rmsdRef <- rmsDist(QCCorr$RefFeats) # 5800887
        rmsdRefClean <- rmsDist(QCCorr$RefFeatsClean) # 3670441
        rmsdRefCorr <- rmsDist(QCCorr$RefFeatsCorr) # 2264153
        rmsdRefFinal <- rmsDist(RefFeatsFinal) # 2209906
        RefRMSD <- data.frame(
            rmsdRef = rmsdRef,
            rmsdRefClean = rmsdRefClean,
            rmsdRefCorr = rmsdRefCorr,
            rmsdRefFinal = rmsdRefFinal
        )
        QCCorr$RefFeatsFinal <- RefFeatsFinal
        QCCorr$RefRMSD <- RefRMSD
    }
    # Write algorithm to remove "removed" variables from clusters
    varClust <- QCCorr$varClust
    N <- length(varClust)
    CVAfter <- nFeatAfter <- numeric(N)
    for (n in seq_len(N)) {
        vars <- varClust[[n]]
        clustVars <- as.matrix(QCFeatsFinal[, finalVars %in% vars])
        CVAfter[n] <- mean(cv(clustVars))
        nFeatAfter[n] <- ncol(clustVars)
    }
    actionInfo <- QCCorr$actionInfo
    actionInfo <- cbind(
        actionInfo[, seq_len(2)],
        nFeatAfter, actionInfo[, 3:5],
        CVAfter
    )
    colnames(actionInfo)[2:3] <- c("nBefore", "nAfter")
    QCCorr$actionInfo <- actionInfo
    QCCorr$QCFeatsFinal <- QCFeatsFinal
    QCCorr$TestFeatsFinal <- TestFeatsFinal
    QCCorr$finalVars <- finalVars
    QCCorr$QCcvs <- QCcvs
    QCFinal <- QCCorr
    message(
        "\nFiltering by QC CV < ", CVlimit, " -> ",
        sum(QCFinal$actionInfo$nAfter), " features out of ",
        sum(QCFinal$actionInfo$nBefore), " kept in the peak table. "
    )
    message("\nPeak table in $TestFeatsFinal, final variables in $finalVars
            and cluster info in $actionInfo.")
    QCFinal
}


#' DC: Perform drift correction
#'
#' Wrapper function for all drift subfunctions 
#' (clust, driftCalc, driftCorr, cleanVar)
#'
#' @param QCObject QC Object (as obtained from `makeQCObject()`)
#' @param modelNames Which MClust geometries to test (see mclust documentation)
#' @param G Which numbers of clusters to test (see mclust documentation)
#' @param BatchObject Batch object (as obtained from `makeBatchObject()`) to be 
#' corrected for drift
#' @param RefObject Optional reference object 
#' (as obtained from `makeBatchObject()`) to validate whether QC correction 
#' improves data quality
#' @param CVlimit QC feature CV limit as final feature inclusion criterion
#' @param report boolean whether to print a pdf report of drift models
#'
#' @return A drift corrected object with QC features CV<limit containing final 
#' peak table
#' @noRd
driftWrap <- function(QCObject,
                        modelNames,
                        G,
                        BatchObject,
                        RefObject = NULL,
                        smoothFunc,
                        spar,
                        CVlimit,
                        report,
                        reportPath) {
    if (is.null(RefObject)) {
        refType <- "none"
    } else {
        refType <- "one"
    }
    driftList <- clust(QCObject$inj,
        QCObject$Feats,
        modelNames = modelNames,
        G = G,
        report = report,
        reportPath = reportPath
    )
    driftList <- driftCalc(driftList,
        smoothFunc = smoothFunc,
        spar = spar,
        report = report,
        reportPath
    )
    driftList <- driftCorr(
        QCDriftCalc = driftList,
        refList = RefObject,
        refType = refType,
        CorrObj = BatchObject,
        report = report,
        reportPath = reportPath
    )
    driftList <- cleanVar(driftList,
        CVlimit = CVlimit,
        report = report,
        reportPath = reportPath
    )
    return(driftList)
}
