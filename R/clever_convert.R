#' Convert \code{"clever_multi"} to \code{"clever"}
#' 
#' @param clev The \code{"clever_multi"} object to convert
#' @return The resulting \code{"clever"} object
#' 
#' @keywords internal
clever_from_multi <- function(clev) {
  class(clev) <- "clever"
  names(clev)[names(clev) == "measures"] <- "measure"
  names(clev)[names(clev) == "outlier_cutoffs"] <- "outlier_cutoff"
  clev$measure_name <- as.character(colnames(clev$measure))
  if (length(clev$measure_name) == 1) {
    # everything except DVARS2
    clev$measure <- as.numeric(clev$measure[,1])
    clev$outlier_cutoff <- as.numeric(clev$outlier_cutoff)
    clev$outlier_flags <- as.logical(clev$outlier_flags[,1])

  }
  clev$ROIs <- NULL
  clev
}

#' Convert \code{"clever"} to \code{"clever_multi"}
#' 
#' Not perfect, since some information has been lost.
#' 
#' @param clev The code{"clever"} object to convert
#' @return The resulting \code{"clever_multi"} object
#' 
#' @keywords internal
clever_to_multi <- function(clev) {
  class(clev) <- "clever_multi"
  names(clev)[names(clev) == "measure"] <- "measures"
  names(clev)[names(clev) == "outlier_cutoff"] <- "outlier_cutoffs"
  if (length(clev$measure_name) == 1) {
    clev$measures <- data.frame(clev$measures)
    clev$outlier_flags <- data.frame(clev$outlier_flags)
    colnames(clev$measures) <- clev$measure_name
    names(clev$outlier_cutoffs) <- clev$measure_name
    colnames(clev$outlier_flags) <- clev$measure_name
  }

  # unresolved: $ROIs, PCA/ICA/PCATF
  
  clev
}
