#' Robust linear model on DCT bases
#' 
#' Fit a linear model regressing an input vector on DCT bases, robustly.
#' 
#' @param x The input vector to regress on DCT bases
#' @param nDCT The number of DCT bases to use. Default: \code{4}
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' 
#' @return The output of \code{robustbase::lmrob}
#' 
#' @keywords internal
rob_trend <- function(x, nDCT=4, lmrob_method="MM") {
  x <- as.vector(x)
  T_ <- length(x)
  
  mat <- data.frame(cbind(1, dct_bases(T_, nDCT)))
  colnames(mat) <- c("x_int", paste0("x_dct", seq(nDCT)))
  mat$y <- x
  
  robustbase::lmrob(y~., mat, method=lmrob_method)
}

#' Variance stabilize a timeseries vector
#' 
#' Variance detrending implemented by the DCT.
#' 
#' @param x The timeseries to variance stabilize
#' @param nDCT The number of DCT bases to use. Default: \code{4}
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' @param rescale After variance stabilizing \code{x}, re-center and re-scale
#'  to the original mean and variance? Default: \code{TRUE}.
#' 
#' @return The variance stabilized timeseries
#' 
#' @export 
var_stabilize <- function(x, nDCT=4, lmrob_method="MM", rescale=TRUE) {
  x_mean <- mean(x); x_var <- var(x)
  x <- as.numeric(scale(x))
  s <- fitted(rob_trend(log(x^2), nDCT, lmrob_method))
  x <- as.numeric(scale(x/s))
  if (rescale) { x <- (x * sqrt(x_var)) + x_mean }
  x
}