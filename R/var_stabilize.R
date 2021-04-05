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

  nDCT <- as.numeric(nDCT)
  stopifnot(nDCT == round(nDCT)); stopifnot(nDCT >= 0)
  if (nDCT == 0) {
    mat <- data.frame(rep(1, T_))
    colnames(mat) <- "x_int"
  } else {
    mat <- data.frame(cbind(1, dct_bases(T_, nDCT)))
    colnames(mat) <- c("x_int", paste0("x_dct", seq(nDCT)))
  }
  mat$y <- x
  
  with(
    set.seed(0),
    robustbase::lmrob(y~., mat, method=lmrob_method, setting="KS2014")
  )
}

#' Variance stabilize a timeseries vector
#'
#' Variance detrending implemented by the DCT.
#'
#' @param x The timeseries to variance stabilize. It should be mean-detrended
#'  already; otherwise, the results will be invalid.
#' @param nDCT The number of DCT bases to use. Default: \code{4}
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' @param rescale After variance stabilizing \code{x}, re-center and re-scale
#'  to the original mean and variance? Default: \code{TRUE}.
#'
#' @return The variance stabilized timeseries
#'
#' @export
#' 
var_stabilize <- function(x, nDCT=2, lmrob_method="MM", rescale=TRUE) {
  if (length(x) < 5) { warning("Timeseries to short to variane stabilize."); return(x) }
  x_mean <- mean(x); x_var <- var(x)
  x <- as.numeric(scale(x))
  s <- as.numeric(rob_trend(log(x^2 + 1), nDCT, lmrob_method)$fitted.values)
  const_mask <- s < 1e-8
  x2 <- x[!const_mask]; s2 <- s[!const_mask]
  s2 <- sqrt(exp(s2) - 1)
  x2 <- x2/s2
  x[!const_mask] <- x2
  if (rescale) { x <- (x * sqrt(x_var)) + x_mean }
  x
}