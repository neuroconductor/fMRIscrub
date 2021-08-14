#' Is this numeric vector constant?
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
check_design_matrix <- function(design, T_=nrow(design)) {
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

#' Cosine bases for the DCT
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

#' Wrapper to common functions for reading NIFTIs
#' 
#' Tries \code{RNifti::readNifti}, then \code{oro.nifti::readNIfTI}. If
#'  neither package is available an error is raised.
#' 
#' @param nifti_fname The file name of the NIFTI.
#' @return The NIFTI
#' @keywords internal
read_nifti <- function(nifti_fname){
  if (requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::readNifti(nifti_fname))
  } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
    return(oro.nifti::readNIfTI(nifti_fname, reorient=FALSE))
  } else {
    stop(
      "Package \"RNifti\" or \"oro.nifti\" needed to read", nifti_fname, 
      ". Please install at least one", call. = FALSE
    )
  }
}

#' Convert to \eqn{T} by \eqn{V} matrix
#' 
#' A \code{"xifti"} is VxT, whereas \code{scrub} accepts a
#'  TxV matrix. This function calls \code{as.matrix} and transposes the data
#'  if it is a \code{"xifti"}.
#' 
#' @param x The object to coerce to a matrix
#' @return x as a matrix.
#' @keywords internal
as.matrix2 <- function(x) {
  if (inherits(x, "xifti")) {
    return( t(as.matrix(x)) )
  } else {
    return( as.matrix(x) )
  }
}