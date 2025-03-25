#' Within-batch signal intensity drift correction
#'
#' Correct drift with cluster-based drift correction. 
#'
#' A basic method for matrix (samples as rows) and SummarizedExperiment is 
#' supported. For grouping variables such as sampleGroups, the basic method 
#' expects vectors, while the SummarizedExperiment method expects the names of 
#' the respective columns. The basic method returns a list with the corrected 
#' peak table and information about the process, whereas the 
#' SummarizedExperiment method assigns the corrected peak table to the object 
#' supplied. 
#'
#' @param peakTable SummarizedExperiment object or matrix
#' @param injections character scalar or numeric, injection order
#' @param sampleGroups character scalar or character, group labels
#' @param QCID character scalar, identifier of QC samples
#' @param RefID character scalar, identifier of external reference 
#' samples for unbiased assessment of drift correction 
#' @param modelNames character, Which mclust geometries to test
#' @param G integer, numbers of clusters to test
#' @param smoothFunc character scalar, choice of regression function: 
#' spline or loess (default: spline)
#' @param spar numeric, smoothing parameter value (defaults to 0.2)
#' @param CVlimit coefficient of variance threshold for filtering 
#' (default = 0.3)
#' @param report boolean, whether to print pdf reports of drift models
#' @param reportPath character scalar, directory path for report
#' @param assay.type character scalar, assay of to be used in case of multiple 
#' assays
#' @param name character scalar, name of the resultant assay in case of multiple
#' assays
#' @param ... parameters with same behavior across methods
#'
#' @return A SummarizedExperiment object with the corrected and filtered matrix 
#' or a list, including the corrected and filtered matrix as well as processing 
#' information:
#' \itemize{
#'   \item actionInfo: see what happened to each cluster
#'   \item testFeatsCorr: to extract drift-corrected data
#'   \item testFeatsFinal: to extract drift-corrected data which pass the 
#'     criterion QC CV < CVlimit
#' }
#'
#' @examples
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
#'     report = FALSE
#' )
#' # Using SummarizedExperiment, more unbiased drift correction using QCs &
#' # external reference samples
#' ## Construct SummarizedExperiment
#' peaks <- SimpleList(t(PTnofill), t(PTfill))
#' sampleData <- meta
#' featureData <- peakInfo(PT = PTnofill, sep = "@", start = 3)
#' rownames(featureData) <- rownames(peaks[[1]])
#' se <- SummarizedExperiment(assays = peaks, colData = sampleData, 
#'                            rowData = featureData)
#' names(assays(se)) <- c("nofill", "fill")
#' 
#' ## Correct drift for single batch
#' se <- se[, colData(se)$batch == "F"]
#' se <- correctDrift(se, 
#'   injections = "inj", sampleGroups = "grp", RefID = "Ref", 
#'   G = seq(5, 35, by = 3), modelNames = c("VVE", "VEE"),
#'   report = FALSE, assay.type = "fill", name = "corrected")
#'
#' @name correctDrift
NULL

.correctDrift <- function(peakTable,
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

#' @export
#' @rdname correctDrift
setGeneric("correctDrift", signature = c("peakTable"),
    function(peakTable, ...) standardGeneric("correctDrift"))

#' @export
#' @rdname correctDrift
setMethod("correctDrift", signature = c("ANY"), .correctDrift)

#' @export
#' @rdname correctDrift
setMethod("correctDrift", signature = c("SummarizedExperiment"),
    function(peakTable, injections, sampleGroups, assay.type, name, ...) {
    from_to <- .get_from_to_names(peakTable, assay.type, name)

    .check_sample_col_present(peakTable, list(injections, sampleGroups))

    # Get corrected peak table and processing metadata
    corrected <- correctDrift(t(assay(peakTable, from_to[[1]])), 
        injections = colData(peakTable)[[injections]],
        sampleGroups = colData(peakTable)[[sampleGroups]], 
        ...)
    corrected_mat <- t(corrected$TestFeatsFinal)
    # Filter peak table by features in corrected peak table, where features
    # with CV < CVlimit are retained
    peakTable <- peakTable[
        which(rownames(peakTable) %in% rownames(corrected_mat)), ]

    # Include corrected peak table in object
    assay(peakTable, from_to[[2]]) <- corrected_mat

    peakTable
})


