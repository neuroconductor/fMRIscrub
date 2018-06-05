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
#' @importFrom gmodels fast.prcomp
#' @examples
#' n_voxels = 1e5
#' n_timepoints = 100
#' x = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' lev = pca_leverage(x)
pca_leverage = function(x,
                        center = TRUE,
                        scale. = FALSE,
                        ...
                        ) {
  res = fast.prcomp(x, center = center,
                    scale. = scale)
}