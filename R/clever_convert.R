#' Convert \code{"clever_multi"} to \code{"clever"}
#'
#' @param clev The \code{"clever_multi"} object to convert
#' @return The resulting \code{"clever"} object
#'
#' @keywords internal
clever_from_multi <- function(clev) {
  class(clev) <- "clever"
  if (ncol(clev$measure) > 1) { 
    warning("The input was not from a `clever` call, since there are more than one leverage measures.")
    return(clev)
  }

  clev$measure <- as.numeric(clev$measure[,1])
  clev$measure_info <- setNames(as.character(clev$measure_info[1,]), colnames(clev$measure_info))
  if ("outlier_cutoff" %in% names(clev)) { clev$outlier_cutoff <- as.numeric(clev$outlier_cutoff) }
  if ("outlier_flag" %in% names(clev)) { clev$outlier_flag <- as.logical(clev$outlier_flag[,1]) }

  # For all projections
  if (grepl("PCA2|ICA2", clev$measure_info["name"])) {
    nComps <- clev$PCA$nPCs_avgvar
  } else {
    nComps <- clev$PCA$nPCs_PESEL
  }
  
  # For PCA
  if (!is.null(clev$PCA$U)) {
    if (nrow(clev$PCA$U) != ncol(clev$PCA$U)) {
      clev$PCA$U <- clev$PCA$U[, seq(nComps), drop=FALSE]
      clev$PCA$D <- clev$PCA$D[seq(nComps), drop=FALSE]
      if ("V" %in% names(clev$PCA)) {
        clev$PCA$V <- clev$PCA$V[, seq(nComps), drop=FALSE]
      }
      if ("highkurt" %in% names(clev$PCA)) {
        clev$PCA$highkurt <- clev$PCA$highkurt[seq(nComps)]
      }
    }
    clev$PCA$nPCs_avgvar <- clev$PCA$nPCs_PESEL <- NULL
  # For PCATF
  } else if (!is.null(clev$PCATF)) {
    if (nrow(clev$PCATF$U) != ncol(clev$PCATF$U)) {
      clev$PCATF$U <- clev$PCATF$U[, seq(nComps), drop=FALSE]
      clev$PCATF$D <- clev$PCATF$D[seq(nComps), drop=FALSE]
      if ("V" %in% names(clev$PCATF)) {
        clev$PCATF$V <- clev$PCATF$V[, seq(nComps), drop=FALSE]
      }
      if ("highkurt" %in% names(clev$PCATF)) {
        clev$PCATF$highkurt <- clev$PCATF$highkurt[seq(nComps)]
      }
    }
    clev$PCA <- NULL
  # For ICA
  } else if (!is.null(clev$ICA)) {
    if (nrow(clev$ICA$M) != ncol(clev$ICA$M)) {
      clev$ICA$S <- clev$ICA$S[, seq(nComps), drop=FALSE]
      clev$ICA$M <- clev$ICA$M[, seq(nComps), drop=FALSE]
      if ("highkurt" %in% names(clev$ICA)) {
        clev$ICA$highkurt <- clev$ICA$highkurt[seq(nComps)]
      }
    }
    clev$PCA <- NULL
  }

  clev
}
