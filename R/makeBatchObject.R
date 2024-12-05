#' Prepare a Batch object for drift correction
#'
#' @param peakTable Peak table (without missing values). Samples in rows; 
#' Features in columns; colnames as mz@rt
#' @param inj Injection sequence number for Batch samples
#' @param QCObject QC Object (to check injection sequence compatibility)
#'
#' @return Batch Object
#'
#' @examples BatchObject <- makeBatchObject(BatchPeakTable,
#'  BatchInjections, QCObject)
#' @noRd
makeBatchObject <- function(peakTable, inj, QCObject) {
    if (length(inj) != nrow(peakTable)) {
        stop("mismatch number of samples in peak table and injection sequence")
    }

    QCInj <- QCObject$inj
    minInj <- min(QCInj)
    maxInj <- max(QCInj)

    if (any(inj < minInj) | any(inj > maxInj)) {
        message("Batch injections outside of QCs. Extrapolation applied.")
    }

    return(list(inj = inj, Feats = peakTable))
}
