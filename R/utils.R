#' Is a numeric vector constant?
#' 
#' @param x The numeric vector
#' @param TOL minimum range of \code{x} to be considered non-constant.
#'  Default: \code{1e-8}
#' 
#' @return Is \code{x} constant? 
#' 
#' @keywords internal
is_constant <- function(x, TOL=1e-8) {
  abs(max(x) - min(x)) < TOL
}

#' Check design matrix
#' 
#' @param design The design matrix
#' 
#' @return The (modified) design matrix
#' 
#' @keywords internal
check_design_matrix <- function(design, T_) {
  class(design) <- "numeric"
  if (identical(design, 1)) { design <- matrix(1, nrow=T_) }
  design <- as.matrix(design)
  stopifnot(nrow(design) == T_)
  # Set constant columns (intercept regressor) to 1, and scale the other columns.
  design_const_mask <- apply(design, 2, is_constant)
  if (any(design_const_mask)) {
    if (any(design_const_mask & abs(design[1,]) < 1e-8)) {
      stop("Constant zero design regressor detected in `design`.")
    }
  }
  design[,design_const_mask] <- 1
  design[,!design_const_mask] <- scale(design[,!design_const_mask])
  design
}

#' Scale data columns robustly
#' 
#' Centers and scales the columns of a matrix robustly for the purpose of 
#'  covariance estimation.
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). If any column MAD is zero, its values become zero
#' and a warning is raised. If all MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#'
#' @return The input matrix with its columns centered and scaled.
#'
#' @importFrom robustbase rowMedians
scale_med <- function(mat){
  TOL <- 1e-8

  # Transpose.
  mat <- t(mat)

  #	Center.
  mat <- mat - c(rowMedians(mat, na.rm=TRUE))

  # Scale.
  mad <- 1.4826 * rowMedians(abs(mat), na.rm=TRUE)
  const_mask <- mad < TOL
  if(any(const_mask)){
    if(all(const_mask)){
    stop("All voxels are zero-variance.\n")
    } else {
      warning(paste0("Warning: ", sum(const_mask),
      " constant voxels (out of ", length(const_mask),
      " ). These will be removed for estimation of the covariance.\n"))
    }
  }
  mad <- mad[!const_mask]
  mat <- mat[!const_mask,]
  mat <- mat/c(mad)

  # Revert transpose.
  mat <- t(mat)

  list(mat=mat, const_mask=const_mask)
}

#' Estimates the parameters of the F distribution of MCD distances.
#'
#' This estimates the parameters c and m required to determine the distribution
#'  of robust MCD distances as derived by Hardin and Rocke (2005), The
#'  Distribution of Robust Distances.
#'
#' @param Q The number of variables in dataset used to compute MCD distances.
#' @param n The total number of observations.
#' @param h The number of observations included in estimation of MCD center and
#'  scale.
#'
#' @return A list containing the estimated F distribution's c, m, and df.
#' @importFrom stats pchisq qchisq
#' @export
fit.F <- function(Q, n, h){
  # Estimate c.
  c <- pchisq(q=qchisq(df=Q, p=h/n), df=Q+2)/(h/n)

  # Estimate asymptotic m.
  alpha <- (n-h)/n
  q_alpha <- qchisq(p=1-alpha, df=Q)
  c_alpha <- (1-alpha)/(pchisq(df=Q+2, q=q_alpha))
  c2 <- -1*pchisq(df=Q+2, q=q_alpha)/2
  c3 <- -1*pchisq(df=Q+4, q=q_alpha)/2
  c4 <- 3*c3
  b1 <- c_alpha*(c3-c4)/(1-alpha)
  b2 <- 0.5 + (c_alpha/(1-alpha))*(c3-q_alpha/Q*(c2 + (1-alpha)/2))
  v1 <- (1-alpha)*b1^2*(alpha*(c_alpha*q_alpha/Q - 1)^2 - 1) -
    2*c3*c_alpha^2*(3*(b1-Q*b2)^2 + (Q+2)*b2*(2*b1-Q*b2))
  v2 <- n*(b1*(b1-Q*b2)*(1-alpha))^2*c_alpha^2
  v <- v1/v2
  m <- 2/(c_alpha^2*v)

  # Corrected m for finite samples.
  m <- m * exp(0.725 - 0.00663*Q - 0.078*log(n))
  df <- c(Q, m-Q+1)

  result <- list(c=c, m=m, df=df)
  return(result)
}

#' Get the cosine bases for DCT
#' 
#' @param T_ Length of timeseries
#' @param n Number of cosine bases
#' 
#' @return Matrix with cosine bases along columns
#' 
#' @export
dct_bases <- function(T_, n){
  b <- matrix(NA, T_, n)
  idx <- (seq(T_)-1)/(T_-1)
  for (ii in seq(n)) { b[,ii] <- cos(idx*pi*ii) }
  b
}

#' Multiple timecourse DCT
#' 
#' Detrends, or estimates the trend, of each column in a matrix using the
#'  discrete cosine transform (DCT).
#'
#' @param x A matrix. (Vectors will be coerced to a matrix with one column.)
#' @param add_mean Add constant column to design matrix to de-mean the data?
#'  Default: \code{TRUE}.
#' @param n_bases Number of cosine bases. Default: \code{4}.
#' @param highpass_freq,TR Not implemented yet; use \code{n_bases}
#' @param residual \code{TRUE} (default) to get the residuals after regressing
#'  the cosine bases (detrend using DCT), or \code{FALSE} to get the estimate
#'  of that regression (trendline using DCT).
#'
#' @return The estimated detrended timeseries, or trend, based on the DCT.
#' 
#' @export
dct_mat <- function(x, add_mean=TRUE, n_bases=4, highpass_freq=NULL, TR=NULL, residual=TRUE){
  x <- as.matrix(x)
  T_ <- nrow(x)

  if (!is.null(highpass_freq)) { stop("Not implemented yet; use `n_bases`") }
  if (!is.null(TR)) { stop("Not implemented yet; use `n_bases`") }
  stopifnot(!is.null(n_bases))

  design <- dct_bases(T_, n_bases)
  # Add constant for mean
  if (add_mean) { design <- cbind(1, design) }
  # Normalize 2-norms to 1, so (X^TX)^(-1) is identity
  design <- design / sqrt((T_+1)/2)
  # X(X^TX)^(-1)X^TB == XX^TB
  xfit <- design %*% t(design) %*% x
  
  if (residual) {
    return(x - xfit)
  } else {
    return(x)
  }
}

#' Wrapper to common functions for reading NIFTIs
#' 
#' @param nifti_fname The file name of the NIFTI.
#' @return The NIFTI.
#' @keywords internal
read_nifti <- function(nifti_fname){
  if (requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::readNifti(nifti_fname))
  } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
    return(oro.nifti::readNIfTI(nifti_fname, reorient=FALSE))
  } else {
    stop("Package \"RNifti\" or \"oro.nifti\" needed to read `X`. Please install at least one", call. = FALSE)
  }
}