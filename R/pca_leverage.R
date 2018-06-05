#' Calculate the Leverage for a Principal Components Analysis
#'
#' @param x data matrix of values, passed to \code{\link{fast.prcomp}}
#'
#' @return A list of components and the leverage
#' @export
#'
#' @importFrom gmodels fast.prcomp
#' @examples
#' x = matrix(rnorm(100*3), ncol = 3)
#' lev = pca_leverage(x)
pca_leverage = function(x) {
  res = fast.prcomp(x)
}