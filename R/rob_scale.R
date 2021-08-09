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

#' Stabilize the center and scale of a timeseries robustly 
#' 
#' Stabilize the center and scale of a timeseries using robust regression of
#'  DCT bases on the first and second moments.
#' 
#' @param x The timeseries to stabilize. 
#' @param center,scale Center and scale? Default: \code{TRUE} for both. If
#'  scaling but not centering, the data must already be centered; otherwise,
#'  the results will be invalid. Can also be the number of DCT bases to use for
#'  robust stabilization of center/scale; \code{TRUE} will use \code{4}.
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' @param rescale After stabilizing \code{x}, re-center and re-scale
#'  to the original mean and variance? Default: \code{TRUE}.
#'
#' @export
#' 
rob_scale <- function(x, center=TRUE, scale=TRUE, lmrob_method="MM", rescale=TRUE) {
  if (length(x) < 5) { warning("Timeseries too short to variance stabilize."); return(x) }
  x_mean <- mean(x); x_var <- var(x)
  x <- as.numeric(scale(x))

  if (isTRUE(center)) {
    center <- 4
  } else if (isFALSE(center)) { 
    center <- 0
  } else {
    center <- as.numeric(center)
    stopifnot(length(center)==1 && center >= 0)
  }

  if (isTRUE(scale)) {
    scale <- 4
  } else if (isFALSE(scale)) { 
    scale <- 0
  } else {
    scale <- as.numeric(scale)
    stopifnot(length(scale)==1 && scale >= 0)
  }

  if (center > 0) {
    m <- as.numeric(rob_trend(x, nDCT=center, lmrob_method)$fitted.values)
    x <- scale(x - m)
  }

  const_mask <- abs(x) < 1e-6

  if (scale > 0) {
    x2 <- ifelse(const_mask, NA, x)
    s <- as.numeric(rob_trend(log(x2^2), nDCT=scale, lmrob_method)$fitted.values)
    s <- sqrt(exp(s))
    if (any(s < 1e-6)) { stop("Error: near-constant variance detected.") } # TEMPORARY
    x[!const_mask] <- x[!const_mask] / s
    x <- scale(x)
  }

  if (rescale) { x <- (x * sqrt(x_var)) + x_mean }
  x
}