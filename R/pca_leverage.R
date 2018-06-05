#' Calculate the Leverage for a Principal Components Analysis
#'
#' @param x data matrix of values, passed to \code{\link{fast.prcomp}}
#' @param center logical value for centering, passed to \code{\link{prcomp}}
#' @param scale. logical value for scaling, passed to \code{\link{prcomp}}
#'
#'
#' @return A list of components and the leverage
#' @export
#'
#' @importFrom gmodels fast.prcomp
#' @examples
#' x = matrix(rnorm(100*3), ncol = 3)
#' lev = pca_leverage(x)
pca_leverage = function(x,
                        center = TRUE,
                        scale. = FALSE
                        ) {
  res = fast.prcomp(x)
}