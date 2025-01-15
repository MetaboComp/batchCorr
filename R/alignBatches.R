#' Between-batch alignment of features
#'
#' Multi-batch alignment of features artificially split between batches.
#' 
#' A basic method with matrices and SummarizedExperiment is supported. For 
#' grouping variables such as batches, the basic method expects vectors, while 
#' the SummarizedExperiment method expects the names of the respective columns 
#' in PeakTabFilled. The basic method returns a list with the aligned peak 
#' table and information about the process, whereas the SummarizedExperiment 
#' method assigns the aligned peak table to the object supplied to 
#' PeakTabFilled. 
#' 
#' @param peakInfo matrix with mz and rt in columns 1:2 (see e.g. ?peakInfo)
#' @param PeakTabNoFill SummarizedExperiment object or matrix before gap-filling
#' @param PeakTabFilled SummarizedExperiment object or matrix without missing 
#' values 
#' @param batches character scalar or factor, batch labels
#' @param sampleGroups character scalar or character, group labels
#' @param selectGroup character scalar, identifier of QC samples
#' @param NAhard numeric scalar, proportion of NAs within batch for feature to 
#' be considered missing
#' @param mzdiff numeric scalar, tolerance for difference in m/z
#' @param rtdiff numeric scalar, tolerance for difference in retention time
#' @param report boolean, whether to export diagnostic plots into your work 
#' directory (default: TRUE)
#' @param reportPath character scalar, directory path for report
#' @param assay.type1 character scalar, assay of PeakTabNoFill to be used in 
#' case of multiple assays
#' @param assay.type2 character scalar, assay of PeakTabFilled to be used in
#' case of multiple assays
#' @param name character scalar, name of the resultant assay in case of 
#' multiple assays
#' @param mz_col character scalar, name of column containing mass information
#' @param rt_col character scalar, name of column for retention time
#' @param ... optional arguments not used
#'
#' @return a SummarizedExperiment object as PeakTabFill with the the aligned 
#' matrix or a list of six:
#' \itemize{
#'   \item PTalign: the aligned peak table
#'   \item boolAveragedAlign: boolean vector of final features where alignment 
#'     has been made using feature averaging
#'   \item PTfill: peaktable without missing values (indata)
#'   \item boolKeep: boolean vector of features kept after alignment
#'   \item boolAveragedFill boolean vector of original features where alignment 
#'     has been made using feature averaging
#'   \item aI: alignIndex object (indata)
#' }
#'
#' @examples
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' }
#' data("ThreeBatchData")
#' # Basic method
#' ## Extract peakinfo (i.e. m/z and rt of features).
#' peakIn <- peakInfo(PT = PTnofill, sep = "@", start = 3)
#' ## Perform multi-batch alignment
#' alignBat <- alignBatches(
#'     peakInfo = peakIn, PeakTabNoFill = PTnofill,
#'     PeakTabFilled = PTfill, batches = meta$batch,
#'     sampleGroups = meta$grp, selectGroup = "QC",
#'     reportPath = "drift_report/"
#' )
#' ## Extract new peak table
#' PT <- alignBat$PTalign
#'
#' # SummarizedExperiment
#' ## Construct SummarizedExperiment
#' peaks <- SimpleList(t(PTnofill), t(PTfill))
#' sampleData <- meta
#' featureData <- peakInfo(PT = PTnofill, sep = "@", start = 3)
#' rownames(featureData) <- rownames(peaks[[1]])
#' se <- SummarizedExperiment(assays = peaks, colData = sampleData, 
#'                            rowData = featureData)
#' names(assays(se)) <- c("nofill", "fill")
#'
#' se <- alignBatches(PeakTabNoFill = se, PeakTabFilled = se,
#'                    batches = "batch", sampleGroups = "grp", report = FALSE,
#'                    assay.type1 = "nofill", assay.type2 = "fill", 
#'                    name = "aligned", rt_col = "rt", mz_col = "mz")
#' \dontshow{
#' setwd(.old_wd)
#' }
#' @name alignBatches

.alignBatches <- function(peakInfo,
                          PeakTabNoFill,
                          PeakTabFilled,
                          batches,
                          sampleGroups,
                          selectGroup = "QC",
                          NAhard = 0.8,
                          mzdiff = 0.002,
                          rtdiff = 15,
                          report = TRUE,
                          reportPath = NULL) {
    if (report & is.null(reportPath)) {
        stop("Argument 'reportPath' is missing")
    }
    if (anyNA(PeakTabFilled)) {
        stop("Missing values in PeakTabFilled")
    }
    if (!selectGroup %in% sampleGroups) {
        stop("Group ", selectGroup, " not in metadata.")
    }
    if (report & is.null(reportPath)) {
        stop("Argument 'reportPath' is missing")
    }
    if (report & !is.null(reportPath)) {
        if (!endsWith(reportPath, "/")) {
            message("Adding slash to path to allow proper folder structure")
            reportPath <- paste0(reportPath, "/")
        }
        if (!file.exists(reportPath)) {
            message("Creating folder ", reportPath)
            dir.create(reportPath, recursive = TRUE)
        }
    }
    bFlag <- batchFlag(
        PTnofill = PeakTabNoFill,
        batch = batches,
        sampleGroup = sampleGroups,
        peakInfo = peakInfo,
        NAhard = NAhard
    )
    aIQ <- alignIndex(
        batchflag = bFlag,
        grpType = selectGroup,
        mzdiff = mzdiff,
        rtdiff = rtdiff,
        report = report,
        reportPath = reportPath
    )
    if (is.null(aIQ)) {
        return(NULL)
    }
    if (report) {
        plotAlign(batchflag = bFlag, alignindex = aIQ, reportPath = reportPath)
    }

    bAlign <- batchAlign(
        batchflag = bFlag,
        alignindex = aIQ,
        peaktable_filled = PeakTabFilled,
        batch = batches
    )
}

setGeneric("alignBatches", signature = c("PeakTabNoFill", "PeakTabFilled"),
    function(PeakTabNoFill, PeakTabFilled, ...) standardGeneric("alignBatches"))

#' @rdname alignBatches
#' @export
setMethod("alignBatches", signature = c("ANY", "ANY"), .alignBatches)

#' @rdname alignBatches
#' @export
setMethod("alignBatches", 
    signature = c(PeakTabNoFill = "SummarizedExperiment", 
                  PeakTabFilled = "SummarizedExperiment"),
    function(PeakTabNoFill, PeakTabFilled, batches, sampleGroups, 
             assay.type1 = NULL, assay.type2 = NULL, name = NULL, 
             mz_col = "mz", rt_col = "rt", ...) {

    from_NoFill <- .get_from_name(PeakTabNoFill, assay.type1)
    # Process name for filled object as it is returned
    from_to_Fill <- .get_from_to_names(PeakTabFilled, assay.type2, name)
  
    .check_sample_col_present(PeakTabNoFill, list(batches, sampleGroups))
    .check_feature_col_present(PeakTabFilled, list(mz_col, rt_col))

    peakIn <- as.matrix(rowData(PeakTabFilled)[, c(mz_col, rt_col)])
  
    aligned <- alignBatches(peakInfo = peakIn, 
        PeakTabNoFill = t(assay(PeakTabNoFill, from_NoFill)),
        PeakTabFilled = t(assay(PeakTabFilled, from_to_Fill[[1]])),
        batches = colData(PeakTabNoFill)[[batches]],
        sampleGroups = colData(PeakTabFilled)[[sampleGroups]],
        ...)
  
    aligned_mat <- t(aligned$PTalign)
  
    # Filter object to accept assay of aligned features
    PeakTabFilled <- PeakTabFilled[which(rownames(PeakTabFilled) %in% 
                                 rownames(aligned_mat)), ]
    assay(PeakTabFilled, from_to_Fill[[2]]) <- aligned_mat
  
    PeakTabFilled
})
