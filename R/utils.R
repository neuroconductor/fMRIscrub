#' Centers and scales a matrix robustly for the purpose of covariance estimation.
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). If any column MAD is zero, its values become zero
#' and a warning is raised. If all MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#'
#' @return The input matrix centered and scaled.
#'
#' @importFrom robustbase rowMedians
scale_med <- function(mat){
  TOL <- 1e-8

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
  N_ <- ncol(mat)

  out <- list(mat=mat, const_mask=const_mask)
  return(out)
}

#' Estimates the trend of \code{ts} using a robust discrete cosine transform.
#'
#' @param ts A numeric vector to detrend.
#' @param robust Should a robust linear model be used? Default FALSE.
#'
#' @return The estimated trend.
#'
#' @importFrom stats mad
#' @importFrom robustbase lmrob
#' @importFrom robustbase lmrob.control
#' @export
est_trend <- function(ts, robust=TRUE){
  EPS <- 1e-8
  if(mad(ts) < EPS){ return(ts) }

  df <- data.frame(
    index=1:length(ts),
    ts=ts
  )

  i_scaled <- 2*(df$index-1)/(length(df$index)-1) - 1 #range on [-1, 1]

  df['p1'] <- cos(2*pi*(i_scaled/4 - .25)) #cosine on [-1/2, 0]*2*pi
  df['p2'] <- cos(2*pi*(i_scaled/2 - .5)) #cosine on [-1, 0]*2*pi
  df['p3'] <- cos(2*pi*(i_scaled*3/4  -.75)) # [-1.5, 0]*2*pi
  df['p4'] <- cos(2*pi*(i_scaled - 1)) # [2, 0]*2*pi

  if(robust){
    control <- lmrob.control(scale.tol=1e-3, refine.tol=1e-2) # increased tol.
    # later: warn.limit.reject=NULL
    trend <- lmrob(ts~p1+p2+p3+p4, df, control=control)$fitted.values
  } else {
    trend <- lm(ts~p1+p2+p3+p4, df)$fitted.values
  }

  return(trend)
}
