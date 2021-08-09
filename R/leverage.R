#' Computes PCA leverage.
#'
#' Computes the leverage of each observation, in the PC score / IC mixing matrix (U/M).
#'  Optionally can identify the outliers.
#'
#' @param Comps The n x Q PC score matrix/IC mixing matrix.
#' @param are_orthogonal Can the columns of \code{Comps} be assumed to be orthogonal
#'  and have a 2-norm of 1? Speeds up the computation.
#' @param median_cutoff The outlier cutoff, in multiples of the median leverage.
#'  Default: \code{NULL} (do not compute outliers).
#' 
#' @return A list with entries \code{"meas"} (the leverage values), 
#'  \code{"cut"} (the leverage cutoff value) and 
#'  \code{"flag"} (logical vector indicating the outliers). If 
#'  \code{!is.null(median_cutoff)}, all entries except \code{"meas"} are omitted.
#'  
#' @importFrom stats median
#' 
#' @export
leverage <- function(Comps, are_orthogonal=FALSE, median_cutoff=NULL){
  if (are_orthogonal) {
    lev <- apply(Comps^2, 1, sum)
  } else {
    #lev <- diag( Comps %*% solve(t(Comps) %*% Comps, t(Comps)) )
    lev <- diag(hat_matrix(Comps))
  }

  out <- list(meas=lev)
  if (!is.null(median_cutoff)){
    out$cut <- median_cutoff * median(lev)
    out$flag <- out$meas > out$cut
  }
  out
}