#' Count each voxel's neighbors with value \code{TRUE}
#'
#' For each voxel in a 3D logical array, count the number of immediate neighbors
#'  with value \code{TRUE}.
#'
#' @param arr The 3D logical array.
#' 
#' @return An array with the same dimensions as \code{arr}. Each voxel value
#'  will be the number of its immediate neighbors (0 to 6) which are \code{TRUE}.
#' 
#' @keywords internal
neighbor_counts <- function(arr){
  stopifnot(length(dim(arr)) == 3)
  stopifnot(all(unique(arr) %in% c(TRUE, FALSE)))
  arrPad <- ciftiTools:::pad_vol(arr, matrix(1, 3, 2), fill=FALSE)
  dPad <- dim(arrPad)
  # Look up, down, left, right, forward, behind (not diagonally)
  arrPad[1:(dPad[1]-2),2:(dPad[2]-1),2:(dPad[3]-1)] +
    arrPad[3:(dPad[1]),2:(dPad[2]-1),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),1:(dPad[2]-2),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),3:(dPad[2]),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),2:(dPad[2]-1),1:(dPad[3]-2)] +
    arrPad[2:(dPad[1]-1),2:(dPad[2]-1),3:(dPad[3])]
}

#' Erode volumetric mask
#'
#' Erode a volumetric mask by a certain number of voxel layers. For each layer,
#'  any in-mask voxel adjacent to at least one out-of-mask voxel is removed
#'  from the mask. 
#'
#' Diagonal voxels are not considered adjacent, i.e. the voxel at (0,0,0) is not
#'  adjacent to the voxel at (1,1,0) or (1,1,1), although it is adjacent to 
#'  (1,0,0).
#'
#' @param vol The volume to erode. Out-of-mask voxels should be indicated by a 
#'  value in \code{out_of_mask_val}.
#' @param n_erosion The number of layers to erode the mask by.
#' @param out_of_mask_val A voxel is not included in the mask if and only if its
#'  value is in this vector. The first value in this vector will be used to replace
#'  the eroded voxels. Default: \code{NA}.
#'
#' @return The eroded \code{vol}. It is the same as \code{vol} but with eroded
#'  voxels replaced with the value \code{out_of_mask_val[1]}.
#'
#' @export
erode_vol <- function(vol, n_erosion=1, out_of_mask_val=NA){
  stopifnot(is.numeric(n_erosion))
  if (n_erosion==0) {return(vol)}
  mask <- !array(vol %in% out_of_mask_val, dim=dim(vol))
  for(ii in 1:n_erosion){
    to_erode <- mask & neighbor_counts(mask) < 6
    mask[to_erode] <- FALSE
    vol[to_erode] <- out_of_mask_val[1]
  }
  vol
}

#' Hat matrix
#'
#' Get the hat matrix of X after adding a column of ones (for the intercept).
#'
#' @param X Numeric matrix. A column of ones will be added for the intercept.
#'
#' @return The hat matrix, a square matrix with the same number of rows/columns
#'  as the number of rows in \code{X}.
#'
#' @keywords internal
hat_mat <- function(X){
  X <- cbind(1, X)
  # The hat matrix is     X * (X^T * X)^(-1) * X^T.
  # This is identical to  X * solve(t(X) %*% X, t(X))
  # Could there be a more efficient way? 
  # https://stackoverflow.com/questions/9071020/compute-projection-hat-matrix-via-qr-factorization-svd-and-cholesky-factoriz/39298028#39298028
  X %*% solve(t(X) %*% X, t(X))
}

#' CompCor
#'
#' The CompCor algorithm (Behzadi et. al., 2007) 10.1016/j.neuroimage.2007.04.042
#' 
#' The data are centered on each voxel timecourse's median.
#'
#' @param vol The \eqn{I x J x K x T} volumetric fMRI timeseries, including both
#'  the data voxels (unless \code{data} is provided) and the noise voxels.
#' @param data_ROI A logical mask for \code{vol} where \code{TRUE} values
#'  indicate data voxels, from which the noise components will be regressed. Up
#'  to one of \code{data_ROI} or \code{data} should be provided.
#' @param data A \eqn{V x T} data matrix from which the noise components
#'  will be regressed. Up to one of \code{data_ROI} or \code{data} should be 
#'  provided.
#' @param noise_ROI A list of logical masks for \code{vol} where \code{TRUE} 
#'  values indicate voxels in a certain "noise ROI" (e.g. the WM or CSF). Noise
#'  components will be obtained separately for each noise ROI.
#' @param noise_nPC A vector the same length as \code{noise_ROI} indicating
#'  the number of PCs from each noise ROI to use. If this vector is a shorter
#'  length its elements will be recycled. The values can also be between 0 and
#'  1, in which case they will represent the minimum proportion of variance 
#'  explained by the PCs used for each noise ROI. The smallest number of PCs
#'  will be used to achieve this proportion of variance explained. Default:
#'  \code{5}.
#' @param noise_erosion The number of layers to erode from each of the 
#'  \code{noise_ROI}. Default: \code{0}.
#'
#' @return A list with entries \code{"data"}, \code{"noise"}, \code{"noise_mask"}
#'  and \code{"noise_var"}.
#'
#'  If \code{data_ROI} or \code{data} was provided, the entry \code{"data"} will be a 
#'  \code{V x T} matrix where each row is a data voxel (if \code{data_ROI} was 
#'  provided, the voxels are in spatial order; if \code{data} was provided, the 
#'  voxels are in the same order) time series with each noise PC regressed from
#'  it. Otherwise, this entry will be \code{NULL}.
#'
#'  The entry \code{"noise"} is a list of \code{T} \eqn{x} \code{noise_nPC} PC 
#'  scores, one for each \code{noise_ROI}.
#'
#' @importFrom robustbase rowMedians
#'
#' @export
CompCor <- function(vol, data_ROI=NULL, data=NULL, noise_ROI, noise_nPC=5, noise_erosion=0){
  stopifnot(length(dim(vol)) == 4)
  T_ <- dim(vol)[4]

  if (!is.null(data_ROI)) {
    stopifnot(dim(data_ROI)==dim(vol)[1:3])
    stopifnot(all(unique(data_ROI) %in% c(TRUE, FALSE)))
    data <- matrix(vol[data_ROI], ncol=T_)
  }
  if (!is.null(data)) {
    stopifnot(is.matrix(data))
    stopifnot(ncol(data)==T_)
  }

  # Get noise components.
  if (is.array(noise_ROI)) { noise_ROI <- list(noise_ROI) }
  stopifnot(is.list(noise_ROI))
  stopifnot(all(noise_nPC > 0))
  N <- length(noise_ROI)
  noise_nPC <- rep(noise_nPC, ceiling(N/length(noise_nPC)))
  noise_erosion <- rep(noise_erosion, ceiling(N/length(noise_erosion)))
  noise_comps <- vector("list", N); names(noise_comps) <- names(noise_ROI)
  noise_var <- vector("list", N); names(noise_var) <- names(noise_ROI)
  for (ii in 1:N) {
    noise_ROI[[ii]] <- erode_vol(noise_ROI[[ii]], noise_erosion[ii], FALSE)
    stopifnot(is.array(noise_ROI[[ii]]) && dim(noise_ROI[[ii]])==dim(vol)[1:3])
    # Vectorize noise ROI.
    noise_ii <- matrix(vol[noise_ROI[[ii]]], ncol=T_)
    noise_ii <- noise_ii - robustbase::rowMedians(noise_ii, na.rm=TRUE)
    # Compute the PC scores.
    if (noise_nPC[ii] >= 1) {
      x <- svd(crossprod(noise_ii), nu=noise_nPC[ii], nv=0)
      noise_comps[[ii]] <- x$u
      noise_var[[ii]] <- ((x$d^2)/sum(x$d^2))[1:noise_nPC[ii]]
    } else {
      x <- svd(crossprod(noise_ii))
      noise_var[[ii]] <- (x$d^2)/sum(x$d^2)
      # Use enough PCs to explain the desired proportion of variance.
      noise_nPC[ii] <- min(length(x$d), sum(cumsum(noise_var[[ii]]) < noise_nPC[ii]) + 1)
      noise_comps[[ii]] <- x$u[,1:noise_nPC[ii]]
      noise_var[[ii]] <- noise_var[[ii]][1:noise_nPC[ii]]
    }
  }

  # Regress.
  if (!is.null(data)) {
    # Project each row of the data (which is why we transpose) on 
    #   the PCs (and the vector of ones for the intercept).
    noise_mat <- do.call(cbind, noise_comps)
    I_minus_H <- diag(nrow(noise_mat)) - hat_mat(noise_mat)
    data <- t(I_minus_H %*% t(data))
  } else {
    data <- NULL
  }

  list(data=data, noise=list(PCs=noise_comps, var=noise_var, ROI=noise_ROI))
}

#' CompCor for HCP data
#'
#' Wrapper to \code{\link{CompCor}} for HCP data. Can perform CompCor on a
#'  data ROI from the same fMRI volume as the one containing the noise ROIs,
#'  or on a separate CIFTI file. 
#'
#' @param vol The \eqn{I x J x K x T} volumetric fMRI timeseries, including both
#'  the data voxels (unless \code{data_cii} is provided) and the noise voxels.
#'  In the HCP, the corresponding file is e.g. 
#'  "Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz". Can also be the file path
#'  to this NIFTI.
#' @param labs The \eqn{I x J x K} labels. Values should be according to this table: 
#'  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT .
#'  In the HCP, the corresponding file is "ROIs/Atlas_wmparc.2.nii.gz". Can also
#'  be the file path to this NIFTI.
#' @param data_labs Numeric vector of labels that make up the data ROI. If 
#'  \code{NULL} (default), use these:
#'
#'  \code{do.call(c, list(cortex=c(1000:2035), subcort=c(8, 10, 11, 12, 13, 16, 17, 18, 26, 28, 47, 49, 50, 51, 52, 53, 54, 58, 60)))}
#'
#'  Alternatively, this can be \code{"cortex"} to use the default labels
#'  corresponding to the cortical grey matter, or \code{"subcort"} to use the
#'  default labels corresponding the the subcortical grey matter.
#'
#'  If \code{data_cii} is provided \code{data_labs} is ignored.
#' @param data_cii \code{"cifti"} object from which the noise ROI
#'  components will be regressed. If \code{data_cii} is provided 
#'  \code{data_labs} is ignored.
#' @param noise_ROI List of numeric vectors, where each entry contains the label
#'  values in \code{labs} corresponding to a noise ROI e.g. the white matter 
#'  (WM) or cerebrospinal fluid (CSF). If \code{NULL} (default), use these:
#'
#'  \code{list(wm_cort = c(3000:4035, 5001, 5002), csf = c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255))}
#'
#'  These default ROIs are based on this forum post: 
#'  https://www.mail-archive.com/hcp-users@humanconnectome.org/msg00931.html
#'
#'  Alternatively, this can be \code{"wm_cort"} to use the default entry
#'  corresponding to cortical white matter, or \code{"csf"} to use the
#'  default entry corresponding to the CSF. 
#' @param noise_nPC,noise_erosion See \code{\link{CompCor}}
#' @param verbose Print occasional updates? Default: \code{TRUE}.
#' @return A list with entries \code{"data"} and \code{"noise"}.

#'  If \code{!is.null(data_labs)}, \code{"data"} will be a \eqn{V x T} data matrix where
#'  each row is a voxel time series from which the noise components have been
#'  regressed (voxels in spatial order, and only voxels in the data ROI are 
#'  included). Otherwise, if \code{!is.null(data_cii)}, it will be the same \code{"cifti"}
#'  object but with the noise components regressed from each brainordinate. 
#'
#'  The entry \code{"noise"} is a list of \code{T} \eqn{x} \code{noise_nPC} PC 
#'  scores, one for each \code{noise_ROI}.
#'
#' @export
CompCor.HCP <- function(vol, labs=NULL, data_labs=NULL, data_cii=NULL, noise_ROI=NULL, noise_nPC=5, noise_erosion=0, verbose=TRUE){
  # Check `vol`.
  if (is.character(vol)) {
    cat("Reading NIFTI.\n")
    if (requireNamespace("RNifti", quietly = TRUE)) {
      vol <- RNifti::readNifti(vol)
    } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
      vol <- oro.nifti::readNIfTI(vol)
    } else {
      stop("Package `RNifti` or `oro.nifti` required to read the NIFTI file. Please install one of them.")
    }
  }
  stopifnot(length(dim(vol))==4)

  # Check `labs`.
  if (is.character(labs)) {
    if (requireNamespace("RNifti", quietly = TRUE)) {
      labs <- RNifti::readNifti(labs)
    } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
      labs <- oro.nifti::readNIfTI(labs)
    } else {
      stop("Package `RNifti` or `oro.nifti` required to read the NIFTI file. Please install one of them.")
    }
  }
  stopifnot(length(dim(labs))==3)
  stopifnot(dim(labs)==dim(vol)[1:3])

  # Check `data_labs` or `data_cii`.
  if (verbose) { cat("Preparing data ROI (to regress from).\n") }
  use_cifti <- !is.null(data_cii)
  if (!use_cifti) {
    data_labs.default <- list(
      cortex = c(1000:2035),
      subcort = c(8, 10, 11, 12, 13, 16, 17, 18, 26, 28, 47, 49, 50, 51, 52, 53, 54, 58, 60)
    )
    if (is.null(data_labs)) { 
      data_labs <- data_labs.default
    } else if (is.character(data_labs)) {
      stopifnot(all(data_labs %in% c("cortex", "subcort")))
      data_labs <- data_labs.default[unique(data_labs)]
    } else if (is.list(data_labs)) {
      data_labs <- data_labs
      stopifnot(is.numeric(do.call(c, data_labs)))
    } else {
      stop("`data_labs` argument is not in a recognized form.")
    }
    data_ROI <- array(labs %in% do.call(c, data_labs), dim=dim(labs))
    data <- list(matrix(vol[data_ROI], ncol=dim(vol)[4]))
  } else {
    if (is.character(data_cii)) { 
      if (requireNamespace("ciftiTools", quietly = TRUE)) {
        data_cii <- ciftiTools::read_cifti(data_cii) 
      } else {
        stop("Package `ciftiTools` required to read the CIFTI file. Please install it from the github repo `mandymejia/ciftiTools`.")
      }
    }
    data <- data_cii$data
    data_ROI <- NULL
  }

  # Check `noise_ROI`.
  if (verbose) { cat("Preparing noise ROIs (to use as regressors).\n") }
  noise_ROI.default <- list(
    wm_cort = c(3000:4035, 5001, 5002), 
    csf = c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255)
  )
  if (is.null(noise_ROI)) {
    noise_ROI <- noise_ROI.default
  } else if (is.character(noise_ROI)) {
    stopifnot(all(noise_ROI %in% c("wm_cort", "csf")))
    noise_ROI <- noise_ROI.default[unique(noise_ROI)]
  } else if (is.list(noise_ROI)) {
    noise_ROI <- noise_ROI
    stopifnot(is.numeric(do.call(c, data_labs)))
  } else {
    stop("`noise_ROI` argument is not in a recognized form.")
  }
  noise_ROI <- lapply(noise_ROI, function(x){array(labs %in% x, dim=dim(labs))})
  x <- CompCor(vol, noise_ROI=noise_ROI, noise_nPC=noise_nPC, noise_erosion=noise_erosion)

  # Project each row of the data (which is why we transpose) on 
  #   the PCs (and the vector of ones for the intercept).
  if (verbose) { cat("Performing CompCor regression.\n") }
  noise_mat <- do.call(cbind, x$noise$PCs)
  I_minus_H <- diag(nrow(noise_mat)) - hat_mat(noise_mat)
  for (ii in 1:length(data)) {
    if (is.null(data[[ii]])) { next }
    data <- t(I_minus_H %*% t(data))
  }

  if (!use_cifti) {
    data <- data[[1]]
  } else {
    data_cii$data <- data
    data <- data_cii
  }

  list(data=data, data_ROI=data_ROI, noise=list(x$noise))
}