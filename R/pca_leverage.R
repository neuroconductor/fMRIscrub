#' Calculate the Leverage for a Principal Components Analysis
#'
#' @param x data matrix of values, passed to \code{\link{fast.prcomp}}
#' @param center logical value for centering, passed to \code{\link{prcomp}}
#' @param scale. logical value for scaling, passed to \code{\link{prcomp}}
#' @param ... additional arguments passed to \code{\link{fast.prcomp}}
#'
#' @return A list of components and the leverage
#' @export
#'
#' @importFrom bootSVD fastSVD
#' @examples
#' n_voxels = 1e5
#' n_timepoints = 100
#' x = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' lev = pca_leverage(x)
pca_leverage = function(
  x,
  center = TRUE,
  scale. = FALSE,
  ...
) {
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  res = bootSVD::fastSVD(x, center_A = FALSE, ...)

  res$sdev = res$d
  res$d = NULL
  res$center = center
  res$scale = scale
  res$rotation = res$v
  res$v = NULL

  class(res) = "prcomp"
  return(res)
}