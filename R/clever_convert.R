#' Convert \code{"clever_multi"} to \code{"clever"}
#'
#' @param clev The \code{"clever_multi"} object to convert
#' @return The resulting \code{"clever"} object
#'
#' @keywords internal
clever_from_multi <- function(clev) {
  class(clev) <- "clever"
  clev$measure_name <- as.character(colnames(clev$measure))

  if (length(clev$measure_name) == 1) {
    # everything except DVARS2
    clev$measure <- as.numeric(clev$measure[,1])
    if ("outlier_cutoff" %in% names(clev)) { clev$outlier_cutoff <- as.numeric(clev$outlier_cutoff) }
    if ("outlier_flags" %in% names(clev)) { clev$outlier_flags <- as.logical(clev$outlier_flags[,1]) }

  }
  clev$ROIs <- NULL

  if ("PCA" %in% names(clev)) {

    # For all projections
    PESEL <- any(grepl("PCA2|ICA2", clev$measure_name))
    nComps <- ifelse(PESEL, clev$PCA$nPCs_PESEL, clev$PCA$nPCs_avgvar)

    # For PCA
    if ("U" %in% names(clev$PCA)) {
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
    } else if ("PCATF" %in% names(clev)) {
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
    } else if ("ICA" %in% names(clev)) {
      if (nrow(clev$ICA$M) != nrow(clev$ICA$M)) {
        clev$ICA$S <- clev$ICA$S[, seq(nComps), drop=FALSE]
        clev$ICA$M <- clev$ICA$M[, seq(nComps), drop=FALSE]
        if ("highkurt" %in% names(clev$ICA)) {
          clev$ICA$highkurt <- clev$ICA$highkurt[seq(nComps)]
        }
      }
      clev$PCA <- NULL
    }
  }

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
  if ("outlier_cutoff" %in% names(clev)) {
    names(clev)[names(clev) == "outlier_cutoff"] <- "outlier_cutoffs"
  }
  if (length(clev$measure_name) == 1) {
    clev$measures <- data.frame(clev$measures)
    colnames(clev$measures) <- clev$measure_name
    if ("outlier_flags" %in% names(clev)) { 
      clev$outlier_flags <- data.frame(clev$outlier_flags)
      colnames(clev$outlier_flags) <- clev$measure_name
    }
    if ("outlier_cutoffs" %in% names(clev)) {
      names(clev$outlier_cutoffs) <- clev$measure_name
    }
  }

  # unresolved: $ROIs, PCA/ICA/PCATF

  clev
}
