#' ThreeBatchData
#'
#' The untargeted LC-MS metabolomics data was collected as per the associated
#' publication. In brief, the data was preprocessed using the \pkg{xcms}
#' package and biological samples were removed.
#'
#' @format A list with three objects encompassing three batches, 48 QC samples,
#' 42 long-term reference samples and 5000 features
#' \describe{
#'   \item{PTnofill}{matrix, peak table with missing values}
#'   \item{PTfill}{matrix, peak table without missing values}
#'   \item{meta}{data.frame, sample metadata including batch, group and
#'   injection order}
#' }
#'
#' @name ThreeBatchData
#' @docType data
#' @author Carl Brunius
#' @references
#' Carl Brunius, Lin Shi, Rikard Landberg
#' Large-scale untargeted LC-MS metabolomics data correction using
#' between-batch feature alignment and cluster-based within-batch signal
#' intensity drift correction.
#' Metabolomics, 12:173. \url{https://doi.org/10.3390/metabo10040135}
NULL

#' @name PTnofill
#' @rdname ThreeBatchData
NULL

#' @name PTfill
#' @rdname ThreeBatchData
NULL

#' @name meta
#' @rdname ThreeBatchData
NULL

###################### SummarizedExperiment helpers ############################
#' @import methods
#' @import SummarizedExperiment
NULL  

.get_from_to_names <- function(object, assay.type, name) {
  # Input behavior (from)
  # If assay.type is not supplied and there is only one assay in the object, 
  # choose the first assay
  if (is.null(assay.type) && length(assays(object)) == 1) {
    assay.type <- 1
  } else if (is.null(assay.type)) {
    stop("When using multiple assays, specify 'assay.type'", call. = FALSE)
  } else if (!assay.type %in% names(assays(object)) & assay.type != 1) {
    stop(assay.type, " was specified but not found in assays", call. = FALSE)
  }
  # Output behavior (to)
  if (is.null(name) && length(assays(object)) == 1) {
    name <- 1
  } else if (is.null(name)) {
    stop("When using multiple assays, specify name of new assay", call. = FALSE)
  } else if (name == assay.type & assay.type != 1) {
    stop("'name' must be different from 'assay.type'", call. = FALSE)
  }
  list(assay.type, name)
}
  
.get_from_name <- function(object, assay.type) {
  # Input behavior (from)
  if (is.null(assay.type) && length(assays(object)) == 1) {
    assay.type <- 1
  } else if (is.null(assay.type)) {
    stop("When using multiple assays, specify assay.type", call. = FALSE)
  } else if (!assay.type %in% names(assays(object)) & assay.type != 1) {
    stop(assay.type, " was specified but not found in assays", call. = FALSE)
  } else {
    assay.type
  }
}

.check_sample_col_present <- function(object, col_names) {
    lapply(col_names, function(col_name) {
        if (is.null(colData(object)[[col_name]])) {
            stop("Column '", col_name, "' not found in colData", call. = FALSE)
        }
    })
}

.check_feature_col_present <- function(object, col_names) {
    lapply(col_names, function(col_name) {
        if (is.null(rowData(object)[[col_name]])) {
            stop("Column '", col_name, "' not found in rowData", call. = FALSE)
        }
    })
}


