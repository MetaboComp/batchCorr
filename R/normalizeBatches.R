#' Between-batch normalisation
#'
#' Batches are feature-wise normalised by Ref samples if passing heuristic 
#' criteria (CV and fold change).
#' Otherwise normalized by population median (still per feature).
#'
#' @param peakTableCorr Multi-batch drift-corrected peak table 
#' (samples in rows; features in columns)
#' @param peakTableOrg Multi-batch original peak table 
#' (samples in rows; features in columns)
#' @param batches Vector (length=nrow(PeakTab)) of batch identifiers
#' @param sampleGroup Vector (length=nrow(PeakTab)) of sample type 
#' (e.g. "sample", "QC", "Ref)
#' @param refGroup Identifier of reference samples in sampleGroups 
#' (defaults to "QC")
#' @param population Identifier of population samples in sampleGroups 
#' (all (default) or any type of samples present in sampleGroups)
#' @param CVlimit  CV criterion to pass for Ref samples per batch 
#' (defaults to 0.3)
#' @param FCLimit Fold-change criterion for intensity 
#' (in relation to average intensity FC between batches; defaults to 5)
#' @param medianZero Strategy for substituting median value for population 
#' normalization when median is zero (-> Inf). Either 'mean' or 'min' 
#' (i.e. lowest non-zero value; default)
#'
#' @return An list object containing:
#' @return peakTable: Normalised multi-batch peak table
#' @return refCorrected: Boolean matrix with info on which batches were 
#' normalised by reference samples; others were normalized by population medians
#'
#' @examples
#' # Note that the example data does not include any biological samples, in
#' # which case population = "sample" does not work
#' data("ThreeBatchData")
#' normData <- normalizeBatches(
#'     peakTableCorr = PTfill, batches = meta$batch,
#'     sampleGroup = meta$grp, refGroup = "Ref",
#'     population = "all"
#' )
#'
#' @importFrom stats median
#'
#' @export
normalizeBatches <- function(peakTableCorr,
                                peakTableOrg,
                                batches,
                                sampleGroup,
                                refGroup = "QC",
                                population = "all",
                                CVlimit = 0.3,
                                FCLimit = 5,
                                medianZero = c("mean", "min")) {
    # Basic sanity checks and setting defaults
    # Setting defaults
    if (missing(peakTableOrg)) {
        peakTableOrg <- peakTableCorr
    }
    if (missing(medianZero)) {
        medianZero <- "min"
    }
    if (population == "all") {
        popSample <- rep(TRUE, nrow(peakTableCorr))
    } else {
        # Checking that user specified population is actually
        # present among the group denotations
        if (!population %in% sampleGroup) {
        stop('population identifier must be present in sampleGroups\nConsider
            setting population="all"')
        } else {
            popSample <- sampleGroup == population
        }
    }
    # If dimensions or colnames of Corr and
    # Org peakTables are not identical throw error
    if (!(identical(dim(peakTableCorr), dim(peakTableOrg)) &
        identical(colnames(peakTableCorr), colnames(peakTableOrg)))) {
        stop("Mismatch between peakTableCorr and peakTableOrg")
    }
    # Checking that refGroup is present in sampleGroup
    if (!refGroup %in% sampleGroup) {
        stop('refGroup identifier needs to be present in sampleGroups.
            \nConsider setting population="all"')
    }
    # Checking that sampleGroup is the same length as batches
    # is the same nrows as peakTableCorr
    if (nrow(peakTableCorr) != length(batches) ||
        nrow(peakTableCorr) != length(sampleGroup)) {
        whichNotCorrectLength <- c("'batches'", "'sampleGroup'")
        stop(paste0(
            paste(
                whichNotCorrectLength[c(
                    nrow(peakTableCorr) != length(batches),
                    nrow(peakTableCorr) != sampleGroup
                )],
                collapse = " & "
            ),
            " not the same length as number of rows in peak table."
        ))
    }

    # Extract info
    uniqBatch <- unique(batches)
    nBatch <- length(uniqBatch)
    nFeat <- ncol(peakTableCorr)

    # Declare/allocate variables
    peakTableNormalized <- peakTableCorr
    # Ref sample CVs per batch (row) and feature (col)
    CVMatrix <- matrix(
        nrow = nBatch,
        ncol = nFeat,
        dimnames = list(uniqBatch, colnames(peakTableCorr))
    )
    # Ref sample mean intensity  per batch (row) and feature (col)
    RefMeanIntensity <- CVMatrix
    # Boolean (flag) for which batches (rows) of which features (columns)
    # are normalized by reference samples
    RefNormMatrix <- matrix(FALSE, nrow = nBatch, ncol = nFeat)
    # Mean intensity ratios between batches (all features)
    MeanIntensityRatios <- matrix(1, nrow = nBatch, ncol = nBatch)


    # Aggregate CV and average intensities for 
    # the reference sample type per batch
    for (b in seq_len(nBatch)) {
        batch <- uniqBatch[b]
        peakTableBatch <- peakTableCorr[batches == batch &
                                            sampleGroup == refGroup, ]
        # Criterion for QC feature CV
        CVMatrix[b, ] <- ifelse(cv(peakTableBatch) <= CVlimit, TRUE, FALSE)
        # Calculate average intensities
        RefMeanIntensity[b, ] <- apply(peakTableBatch, 2, mean)
    }

    # Calculate average intensity ratios between batches
    for (b in seq_len(nBatch - 1)) {
        for (bb in (b + 1):nBatch) {
            MeanIntensityRatios[bb, b] <- 
                mean(RefMeanIntensity[bb, ]) / mean(RefMeanIntensity[b, ])
            
            MeanIntensityRatios[b, bb] <- 1 / MeanIntensityRatios[bb, b]
        }
    }


    # Perform normalisation per feature
    for (feat in seq_len(nFeat)) {
        # Calculate feature-wise average intensity ratios between batches
        featureIntensityRatios <- matrix(1, nrow = nBatch, ncol = nBatch)
        for (b in seq_len(nBatch - 1)) {
            for (bb in (b + 1):nBatch) {
                featureIntensityRatios[bb, b] <- 
                    RefMeanIntensity[bb, feat] / RefMeanIntensity[b, feat]
                featureIntensityRatios[b, bb] <- 
                    1 / featureIntensityRatios[bb, b]
            }
        }

        # Identify candidates for ref normalization
        # Criterion for intensity ratio
        candidates <- abs(log(featureIntensityRatios / MeanIntensityRatios)) <= 
            log(FCLimit)
        
        # Convert missing values -> FALSE (i.e. not a candidate)
        candidates[is.na(candidates)] <- FALSE
        # Criterion for CV < limit
        for (b in seq_len(nBatch)) {
            if (CVMatrix[b, feat] == FALSE | is.na(CVMatrix[b, feat])) {
                candidates[, b] <- candidates[b, ] <- FALSE
            }
        }

        # Perform default normalization by reference samples 
        # within the selected sample populations
        
        # Find reference batch
        refBatch <- min(which(colSums(candidates) == max(colSums(candidates))))
        # Extract flags for reference sample correction
        refCorrFlags <- candidates[, refBatch] 
        # Find which batches to normalise to the reference
        refCorrIndex <- which(refCorrFlags) 
        refCorrIndex <- refCorrIndex[refCorrIndex != refBatch]

        refIntensity <- RefMeanIntensity[refBatch, feat] # Ref batch intensity
        for (b in refCorrIndex) { # Correct batches to reference batch intensity

            # Calculate correction factor for "population" samples and use for
            # normalizaiton
            correctionFactor <- refIntensity / RefMeanIntensity[b, feat]
            peakTableNormalized[which(batches == uniqBatch[b]), feat] <- 
                peakTableCorr[which(batches ==
                                        uniqBatch[b]), feat] * correctionFactor
        }

        # Store "flag" of whether the feature was normalized by ref samples
        RefNormMatrix[, feat] <- refCorrFlags

        # Perform population (median) normalisation for any other batches

        # Check if any batches were not normalized by the reference samples
        if (length(refCorrIndex) + 1 != nBatch) {
            # All samples of "population" corrected by reference samples
            refCorrected <- 
                peakTableNormalized[batches %in%
                                        uniqBatch[c(refBatch, refCorrIndex)] &
                                            popSample, feat]
            refCorrMedian <- median(refCorrected) # Extract their median value

            # Which batches to correct by population median instead
            WhichPOPCorr <- which(!refCorrFlags)

            # Correct those by population median approach
            for (n in WhichPOPCorr) {
                populationMedian <- 
                    median(peakTableOrg[which(batches == uniqBatch[n] &
                                                  popSample), feat])

                # "Fix" for if population median == 0, which will
                # cause division by zero
                if (populationMedian == 0) {
                    if (medianZero == "min") {
                        populationMedian <- 
                                peakTableOrg[which(batches == uniqBatch[n] &
                                                       popSample), feat]
                        populationMedian <- 
                            ifelse(sum(populationMedian != 0) > 0,
                            min(populationMedian[populationMedian != 0]),
                            0
                        )
                    } else if (medianZero == "mean") {
                        populationMedian <- 
                            mean(peakTableOrg[which(batches == uniqBatch[n] &
                                                        popSample), feat])
                    } else {
                        stop("Other options not included at present.")
                    }
                }

                # Calculate corr factor for "population" samples 
                # and use for normalization
                correctionFactor <- ifelse(populationMedian != 0,
                                        refCorrMedian / populationMedian,
                                        1)
                peakTableNormalized[which(batches == uniqBatch[n]), feat] <- 
                    peakTableOrg[which(batches == uniqBatch[n]), feat] *
                    correctionFactor
            }
        }
    }

    return(list(
        peakTable = peakTableNormalized,
        refCorrected = RefNormMatrix
    ))
}
