#' Within-batch signal intensity drift correction
#'
#' @param peakTable Batch peak table
#' @param injections Injection number in sequence
#' @param sampleGroups Vector (length=nrow(PeakTab)) of sample type 
#' (e.g. "sample", "QC", "Ref)
#' @param QCID QC identifier in sampleGroups (defaults to "QC")
#' @param RefID Optional identifier in sampleGroups of external reference 
#' samples for unbiased assessment of drift correction 
#' (see batchCorr tutorial at Gitlab page)
#' @param modelNames Which MClust geometries to test (see mclust documentation)
#' @param G Which numbers of clusters to test (see mclust documentation)
#' @param smoothFunc choice of regression function: spline or loess 
#' (defaults to spline)
#' @param spar smoothing parameter for spline or loess (defaults to 0.2)
#' @param CVlimit QC feature CV limit as final feature inclusion criterion
#' @param report Boolean whether to print pdf reports of drift models
#' @param reportPath directory path for report
#'
#' @return A driftCorrection object
#' @return $actionInfo (to see what happened to each cluster)
#' @return $testFeatsCorr (to extract drift-corrected data) and
#' @return $testFeatsFinal (to extract drift-corrected data which 
#' pass the criterion that QC CV < CVlimit.
#' @export
#'
#' @examples
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' }
#' data("ThreeBatchData")
#' set.seed(2024)
#' # Get batches
#' batchB <- getBatch(
#'     peakTable = PTfill, meta = meta,
#'     batch = meta$batch, select = "B"
#' )
#' batchF <- getBatch(
#'     peakTable = PTfill, meta = meta,
#'     batch = meta$batch, select = "F"
#' )
#' # Drift correction using QCs
#' BCorr <- correctDrift(
#'     peakTable = batchB$peakTable,
#'     injections = batchB$meta$inj,
#'     sampleGroups = batchB$meta$grp, QCID = "QC",
#'     G = seq(5, 35, by = 3), modelNames = c("VVE", "VEE"),
#'     reportPath = "drift_report/"
#' )
#' # More unbiased drift correction using QCs & external reference samples
#' FCorr <- correctDrift(
#'     peakTable = batchF$peakTable,
#'     injections = batchF$meta$inj,
#'     sampleGroups = batchF$meta$grp, QCID = "QC",
#'     RefID = "Ref", G = seq(5, 35, by = 3),
#'     modelNames = c("VVE", "VEE"),
#'     reportPath = "drift_report/"
#' )
#' # Merge batches for batch normalization, for example
#' mergedData <- mergeBatches(list(BCorr, FCorr))
#' \dontshow{
#' setwd(.old_wd)
#' }
correctDrift <- function(peakTable,
                            injections,
                            sampleGroups,
                            QCID = "QC",
                            RefID = "none",
                            modelNames = c("VVV",
                                            "VVE",
                                            "VEV", 
                                            "VEE", 
                                            "VEI",
                                            "VVI", 
                                            "VII"),
                            G = seq(5, 35, by = 10),
                            smoothFunc = "spline",
                            spar = 0.2,
                            CVlimit = 0.3,
                            report = TRUE,
                            reportPath = NULL) {
    # Some basic sanity check
    if (nrow(peakTable) != length(injections)) {
        stop("nrow(peakTable) not equal to length(injections)")
    }
    if (is.null(colnames(peakTable))) {
        stop("All features/variables need to have unique names")
    }
    if (length(sampleGroups) != length(injections)) {
        stop("length(sampleGroups) not equal to length(injections)")
    }
    if (!identical(sort(injections), injections)) {
        stop("injection sequence is not in order\nPlease resort peakTable,
            injections and sampleGroups accordingly")
    }
    if (report & is.null(reportPath)) {
        stop("Argument 'reportPath' is missing")
    }
    if (report & !is.null(reportPath)) {
        if (!endsWith(reportPath, "/")) {
        message("Adding a slash to file path to allow proper folder structure")
            reportPath <- paste0(reportPath, "/")
        }
        if (!file.exists(reportPath)) {
            message("Creating folder ", reportPath)
            dir.create(reportPath, recursive = TRUE)
        }
    }
    meta <- data.frame(injections, sampleGroups)
    # Extract QC info
    batchQC <- .getGroup(
        peakTable = peakTable,
        meta = meta,
        sampleGroup = sampleGroups,
        select = QCID
    )
    # Prepare QC object for drift correction
    QCObject <- makeQCObject(
        peakTable = batchQC$peakTable,
        inj = batchQC$meta$inj
    )
    # Prepare batch data
    BatchObject <- makeBatchObject(
        peakTable = peakTable,
        inj = injections,
        QCObject = QCObject
    )
    # Prepare external reference data
    if (RefID != "none") {
        # Prepare ref object for drift correction
        batchRef <- .getGroup(
            peakTable = peakTable,
            meta = meta,
            sampleGroup = sampleGroups,
            select = RefID
        )
        RefObject <- makeBatchObject(
            peakTable = batchRef$peakTable,
            inj = batchRef$meta$inj,
            QCObject = QCObject
        )
        # Perform drift correction
        Corr <- driftWrap(
            QCObject = QCObject,
            BatchObject = BatchObject,
            RefObject = RefObject,
            modelNames = modelNames,
            G = G,
            smoothFunc = smoothFunc,
            spar = spar,
            CVlimit = CVlimit,
            report = report,
            reportPath = reportPath
        )
    } else {
        # Perform drift correction
        Corr <- driftWrap(
            QCObject = QCObject,
            BatchObject = BatchObject,
            modelNames = modelNames,
            G = G,
            smoothFunc = smoothFunc,
            spar = spar,
            CVlimit = CVlimit,
            report = report,
            reportPath = reportPath
        )
    }
    return(Corr)
}
