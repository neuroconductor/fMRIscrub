#' Calculate the leverage images
#'
#' @param clev A "clever" object.
#' @param timepoints The times for which to compute leverage images (rows of U).
#' @param const_mask Mask that is TRUE where voxels were removed.
#'
#' @return A list of three: the mean leverage image for each outlier meeting
#'  the thresold, the top leverage image for each outlier, and the indices of
#'  the top leverage images.
#'
# TEMPORARY: internal
#' @keywords internal
get_leverage_images <- function(clev, timepoints=NULL, const_mask=NULL){

  if (is.null(timepoints)) {
    timepoints <- which(clev$outlier_flags)
    if (!(length(timepoints) > 0)) {
      stop("`timepoints=NULL` will get leverage images for outliers, but no outliers detected.")
    }
  } else {
    stopifnot(length(timepoints) > 0)
  }

  if ("PCA" %in% names(clev)) {
    U <- clev$PCA$U
    if (!("V" %in% names(clev$PCA))) { 
      stop("No directions. Run clever again with `solve_dirs=TRUE`.") 
    }
    V <- clev$PCA$V
  } else if ("PCATF" %in% names(clev)) {
    U <- clev$PCATF$U
    if (!("V" %in% names(clev$PCA))) { 
      stop("No directions. Run clever again with `solve_dirs=TRUE`.")
    }
    V <- clev$PCATF$V
  } else if ("ICA" %in% names(clev)) {
    U <- clev$ICA$M
    V <- clev$ICA$S
  }

  stopifnot(all(timepoints %in% seq(nrow(U))))

  if(is.null(const_mask)){ const_mask = rep(FALSE, nrow(V)) }
  N_ <- length(const_mask)
  n_imgs <- length(timepoints)

  lev_imgs <- list(mean=matrix(NA, nrow=n_imgs, ncol=N_))
  lev_imgs$mean[,!const_mask] <- U[timepoints,] %*% t(V)

  lev_imgs$top <- matrix(NA, nrow=n_imgs, ncol=N_)
  lev_imgs$top_dir <- vector(mode="numeric", length=n_imgs)
  for(i in 1:n_imgs){
    idx <- timepoints[i]
    lev_imgs$top_dir[i] <- which.max(U[idx,])[1]
    lev_imgs$top[i,!const_mask] <- V[,lev_imgs$top_dir[i]] #Tie: use PC w/ more var.
  }

  row.names(lev_imgs$mean) <- timepoints
  row.names(lev_imgs$top) <- timepoints
  names(lev_imgs$top_dir) <- timepoints

  lev_imgs
}

#' Applies a 2D/3D mask to a matrix to get a 3D/4D volume time series.
#' @param mat A matrix whose rows are observations at different times, and
#'  columns are pixels/voxels.
#' @param mask A corresponding binary mask, with 1's representing regions
#'  within the area of interest and 0's representing regions to mask out.
#' @param out_of_mask_value Fill value for out-of-mask voxels. Default: \code{NA}.
#' @param sliced_dim If the mask is 2D, which dimension does it represent?
#'  Will default to the 3rd dimension (axial).
#'
#' @return A 4D array representing the volume time series. Time is on the 4th
#'  dimension.
#'
#' @export
Matrix_to_VolumeTimeSeries <- function(mat, mask, out_of_mask_value=NA, sliced_dim = NA){
  in_mask <- mask > 0
  T_ <- nrow(mat)

  if(length(dim(mask)) == 3){
    dims <- c(dim(mask), T_)
  } else if(length(dim(mask)) == 2) {
    if(is.na(sliced_dim)){ sliced_dim=3 } #default to 3rd dim (axial)
    dims <- switch(sliced_dim,
                   c(1, dim(mask), T_),
                   c(dim(mask)[1], 1, dim(mask)[2], T_),
                   c(dim(mask), 1, T_)
    )
  } else {
    stop("Not Implemented: mask must be 2D or 3D.")
  }

  vts <- array(out_of_mask_value, dim=dims)
  for(i in 1:T_){
    vts[,,,i][in_mask] <- mat[i,]
  }

  return(vts)
}

#'  \item{lev_images}{
#'    \describe{
#'      \item{mean}{The average of the PC directions, weighted by the unscaled
#'        PC scores at each outlying time point (U[i,] * V^T). Row names are
#'        the corresponding time points.}
#'      \item{top}{The PC direction with the highest PC score at each outlying
#'        time point. Row names are the corresponding time points.}
#'      \item{top_dir}{The index of the PC direction with the highest PC score
#'        at each outlying time point. Named by timepoint.}
#'    }
#'  }