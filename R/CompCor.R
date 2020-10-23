#' Count each voxel's neighbors with value \code{TRUE}
#'
#' For each voxel in a 3D logical array, count the number of immediate neighbors
#'  with value \code{TRUE}.
#'
#' @param arr The 3D logical array.
#' @param pad Pad value for edge.
#' 
#' @return An array with the same dimensions as \code{arr}. Each voxel value
#'  will be the number of its immediate neighbors (0 to 6) which are \code{TRUE}.
#' 
#' @keywords internal
neighbor_counts <- function(arr, pad=FALSE){
  stopifnot(length(dim(arr)) == 3)
  stopifnot(all(unique(arr) %in% c(TRUE, FALSE)))
  arrPad <- ciftiTools:::pad_vol(arr, matrix(1, 3, 2), fill=pad)
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
    to_erode <- mask & neighbor_counts(mask, pad=TRUE) < 6
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

#' Format data for clever and CompCor
#'
#' @inheritParams data_clever_CompCor_Params
#' @param noise_nPC The number of principal components to compute for each noise
#'  ROI. Alternatively, values between 0 and 1, in which case they will 
#'  represent the minimum proportion of variance explained by the PCs used for
#'  each noise ROI. The smallest number of PCs will be used to achieve this 
#'  proportion of variance explained. 
#' 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its
#'  elements will be recycled. Default: \code{5} (compute the top 5 PCs for 
#'  each noise ROI).
#' @param noise_erosion The number of voxel layers to erode the noise ROIs by. 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its 
#'  elements will be recycled. Default: \code{NULL}, which will use a value of
#'  0 (do not erode the noise ROIs).
#'
#' @return A list with components "X", "X_noise", "ROI_data", and "ROI_noise"
#'
#' @keywords internal
format_data_clever_CompCor <- function(X, ROI_data=NULL, ROI_noise=NULL, noise_nPC=5, noise_erosion=NULL){

  # if X is a file, read it.
  if (is.character(X)) {
    if (endsWith(X, ".dtseries.nii") | endsWith(X, ".dscalar.nii")) {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read `X`. Please install it", call. = FALSE)
      }
      X <- read_cifti(X, brainstructures="all")
    } else {
      X <- read_nifti(X)
    }
  } 

  # X must be a matrix, array, or "xifti"
  if (is.matrix(X)) {
    T_ <- nrow(X); V_ <- ncol(X)
    if (T_ > V_) {
      warning(
        "Data matrix has more rows than columns. Check that observations\
        are in rows and variables are in columns."
      )
    }
    X_type <- "vector"
  } else if (is.array(X)) {
    if (length(dim(X))==3) { X <- array(X, dim=c(dim(X)[1:2], 1, dim(X)[3])) }
    stopifnot(length(dim(X))==4)
    T_ <- dim(X)[4]
    X_type <- "volume"
  } else if (inherits(X, "xifti")) {
    xifti_meta <- X$meta
    X <- t(do.call(rbind, X$data))
    T_ <- nrow(X); V_ <- ncol(X)
    X_type <- "xifti"
  } else {
    stop("`X` must be a matrix, array, NIFTI, path to a NIFTI, CIFTI, or path to a CIFTI.")
  }

  if (T_ < 2) { stop("There are less than two timepoints.") }

  # ROI_data
  ROI_data_was_null <- is.null(ROI_data)
  if (X_type == "vector") {
    if (ROI_data_was_null) { ROI_data <- rep(TRUE, V_) }
    ROI_data <- as.vector(ROI_data)
    if (length(ROI_data) != V_) { 
      stop("The `ROI_data` must be a logical vector with the same length as columns in the data.") 
    }
    ROI_data <- as.logical(ROI_data)
  } else if (X_type == "volume") {
    if (ROI_data_was_null) { ROI_data <- c(0, NA, NaN) }
    if (is.character(ROI_data)) {
      ROI_data <- read_nifti(ROI_data)
    }
    if (is.vector(ROI_data)) {
      if (length(ROI_data) != length(unique(ROI_data))) {
        stop(
          "`X` is a volume, and `ROI_data` is a vector, so `ROI_data` should be\
          values which out-of-mask voxels take on. But, the values of `ROI_data`\
          are not unique."
        )
      }
      ROI_data <- apply(array(X %in% ROI_data, dim=dim(X)), 1:3, function(x){!all(x)})
    } else if (is.array(ROI_data)) {
      if (length(dim(ROI_data))==2) { ROI_data <- array(ROI_data, dim=c(dim(ROI_data), 1)) }
      if (all(dim(ROI_data) != dim(X)[1:3])) { 
        stop("The `ROI_data` must have the same dimensions as the first three dimensions of `X`.")
      }
      ROI_data[,,] <- as.logical(ROI_data)
    }
  } else if (X_type == "xifti") {
    ROI_data <- rep(TRUE, V_)
  } else { stop("Internal error: unrecognized `X_type`") }
  stopifnot(sum(ROI_data) > 0)

  # ROI_noise
  if (!is.null(ROI_noise)) {
    if (!is.list(ROI_noise)) { ROI_noise <- list(Noise1=ROI_noise) }
    if (is.null(names(ROI_noise))) { names(ROI_noise) <- paste0("Noise", 1:length(ROI_noise)) }
    if (length(names(ROI_noise)) != length(unique(names(ROI_noise)))) {
      stop("The `ROI_noise` names must be unique.")
    }
    stopifnot(!any(names(ROI_noise) == "data"))

    # noise_nPC
    noise_nPC <- as.list(noise_nPC)
    if (is.null(names(noise_nPC))) {
      noise_nPC <- noise_nPC[rep(1:length(noise_nPC), length(ROI_noise))[1:length(ROI_noise)]]
      names(noise_nPC) <- names(ROI_noise)
    }     
    else {
      if (all(sorted(names(noise_nPC)) != sorted(names(ROI_noise)))) {
        stop("The names of `noise_nPC` do not match those of `ROI_noise`.")
      }
    }  
    noise_nPC <- noise_nPC[names(ROI_noise)]

    # noise_erosion
    if (X_type == "volume") {
      if (is.null(noise_erosion)) { 
        noise_erosion = 0
      } else {
        if (!all(noise_erosion==0) & !any(sapply(ROI_noise, is.array))) {
          warning("`noise_erosion` was provided, but there are no array/NIFTI noise ROIs to erode.")
        }
      }
      noise_erosion <- as.list(noise_erosion)
      if (is.null(names(noise_erosion))) {
        noise_erosion <- noise_erosion[rep(1:length(noise_erosion), length(ROI_noise))[1:length(ROI_noise)]]
        names(noise_erosion) <- names(ROI_noise)
      } else {
        if (all(sorted(names(noise_erosion)) != sorted(names(ROI_noise)))) {
          stop("The names of `noise_erosion` do not match those of `ROI_noise`.")
        }
      }
      noise_erosion <- noise_erosion[names(ROI_noise)]
      
    } else {
      if (!is.null(noise_erosion)) { 
        warning(
          "Erosion requires volumetric data, but the data is not volumetric.\
          No erosion will happen."
        ) 
      }
    }

    X_noise <- vector("list", length(ROI_noise)); names(X_noise) <- names(ROI_noise)
    for (ii in 1:length(ROI_noise)) {
      if (is.null(ROI_noise[[ii]])) { ROI_noise[[ii]] <- NULL; next }
      if (X_type == "vector") {
        if (is.vector(ROI_noise[[ii]])) {
          stopifnot(length(ROI_noise[[ii]]) == V_)
          ROI_noise[[ii]] <- as.logical(ROI_noise[[ii]])
          X_noise[[ii]] <- X[,ROI_noise[[ii]]]
        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
        } else {
          stop(
            "Each entry in `ROI_noise` must be a logical vector, or matrix\
            with the same number of rows as timepoints in `X`."
          )
        }
      } else if (X_type == "volume") {
        if (is.character(ROI_noise[[ii]])) {
          if (!file.exists(ROI_noise[[ii]])) { 
            stop(paste(
              "The `ROI_noise` entry", ROI_noise[[ii]], "is not an existing file."
            )) 
          }
          ROI_noise[[ii]] <- read_nifti(ROI_noise[[ii]])
        }
        if (is.matrix(ROI_noise[[ii]])) { 
          ROI_noise[[ii]] <- array(ROI_noise[[ii]], dim=c(dim(ROI_noise[[ii]], 1)))
        }
        if (is.array(ROI_noise[[ii]])) {
          stopifnot(all(dim(ROI_noise[[ii]]) == dim(X)[1:3]))
          ROI_noise[[ii]][,,] <- as.logical(ROI_noise[[ii]]) * 1
          ROI_noise[[ii]] <- erode_vol(ROI_noise[[ii]], noise_erosion[[ii]], c(-1, 0, NA))
          X_noise[[ii]] <- matrix(X[ROI_noise[[ii]] > 0], ncol=T_)

        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
        } else {
          stop(
            "Each entry in `ROI_noise` must be a logical array, or matrix\
            with the same number of rows as timepoints in `X`."
          ) 
        }
      } else if (X_type == "xifti") {
        stopifnot(is.matrix(ROI_noise[[ii]]))
        stopifnot(nrow(ROI_noise[[ii]]) == T_)
        X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
      } else { stop("Internal error: unrecognized `X_type`") }
      if (ncol(X_noise[[ii]]) == 0) { 
        warning(paste("The noise ROI", names(ROI_noise)[ii], "is empty."))
      }
    }

    if (!all(sapply(ROI_noise, is.null))) {
      # check that ROI are mutually exclusive
      all_ROI_noises <- apply(do.call(rbind, ROI_noise) > 0, 2, sum)
      if (!all(all_ROI_noises < 2)) {
        stop("The noise ROIs must all be mutually exclusive.")
      }
      all_ROI_noises <- all_ROI_noises > 0
      if (ROI_data_was_null) { 
        if (X_type == "volume") {
          ROI_data[,,][all_ROI_noises] <- FALSE
        } else {
          ROI_data[all_ROI_noises] <- FALSE
        }
      } else {
        if (any(all_ROI_noises & as.vector(ROI_data))) {
          stop("The noise ROIs must not overlap with the data ROI.")
        }
      }
    }

    if (X_type == "volume") {
      X <- t(matrix(X[!all_ROI_noises], ncol=dim(X)[4]))
    } else {
      X <- X[,!all_ROI_noises]
    }

  } else {
    X_noise <- NULL
  }

  # make sure nPCs greater than the rank of the noise ROIs!
  # similarly for data
  # ...

  list(
    X=X, X_noise=X_noise, 
    ROI_data=ROI_data, ROI_noise=ROI_noise, 
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )
}

#' CompCor: get noise components
#'
#' @param X_noise The noise ROIs data
#' @param noise_nPC Number of PCs to obtain for each noise ROI
#'
#' @return A list with components X, X_noise, ROI_data, ROI_noise, noise_nPC,
#'  noise_erosion, noise_comps, and noise_var.
#' 
#' @keywords internal
CompCor.noise_comps <- function(X_noise, noise_nPC){

  noise_comps <- vector("list", N); names(noise_comps) <- names(X_noise)
  noise_var <- vector("list", N); names(noise_var) <- names(X_noise)

  for (ii in 1:N) {
    if (ncol(X_noise[[ii]]) == 0) { next }
    X_noise[[ii]] <- t(t(X_noise[[ii]]) - robustbase::rowMedians(t(X_noise[[ii]]), na.rm=TRUE))
    # Compute the PC scores.
    if (noise_nPC[ii] >= 1) {
      x <- svd(tcrossprod(X_noise[[ii]]), nu=noise_nPC[ii], nv=0)
      noise_comps[[ii]] <- x$u
      noise_var[[ii]] <- ((x$d^2)/sum(x$d^2))[1:noise_nPC[ii]]
    } else {
      x <- svd(tcrossprod(X_noise[[ii]]))
      noise_var[[ii]] <- (x$d^2)/sum(x$d^2)
      # Use enough PCs to explain the desired proportion of variance.
      noise_nPC[ii] <- min(length(x$d), sum(cumsum(noise_var[[ii]]) < noise_nPC[ii]) + 1)
      noise_comps[[ii]] <- x$u[,1:noise_nPC[ii]]
      noise_var[[ii]] <- noise_var[[ii]][1:noise_nPC[ii]]
    }
  }

  list(noise_comps=noise_comps, noise_var=noise_var)
}

#' CompCor: regress
#'
#' @param X The data
#' @param noise_comps The noise components
#'
#' @return The data with the noise components regressed from it.
#' 
#' @keywords internal
CompCor.regress <- function(X, noise_comps){
  # Project each row of the data (which is why we transpose) on 
  #   the PCs (and the vector of ones for the intercept).
  noise_mat <- do.call(cbind, noise_comps)
  I_minus_H <- diag(nrow(noise_mat)) - hat_mat(noise_mat)
  t(I_minus_H %*% t(X))
}

#' CompCor
#'
#' The CompCor algorithm (Behzadi et. al., 2007) 10.1016/j.neuroimage.2007.04.042
#' 
#' The data are centered on each voxel timecourse's median.
#'
#' @inheritParams data_clever_CompCor_Params
#' @inheritParams noise_Params
#'
#' @return A list with entries \code{"data"} and \code{"noise"}
#'
#'  The entry \code{"data"} will be a \code{V x T} matrix where each row is a 
#'  data voxel (if it was originally an array, the voxels will be in spatial
#'  order) time series with each noise PC regressed from it. Otherwise, this 
#'  entry will be \code{NULL}.
#'
#'  The entry \code{"noise"} is a list of \code{T} \eqn{x} \code{noise_nPC} PC 
#'  scores, one for each \code{ROI_noise}.
#'
#' @importFrom robustbase rowMedians
#'
#' @export
CompCor <- function(X, ROI_data=NULL, ROI_noise=NULL, noise_nPC=5, noise_erosion=NULL){

  out1 <- format_data_clever_CompCor(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise, noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )

  out2 <- CompCor.noise_comps(
    X_noise=out1$X_noise, noise_nPC=out1$noise_nPC
  )

  out1$X <- CompCor.regress(out1$X, out2$noise_comps)

  list(
    X=out1$X, 
    noise=list(PCs=out2$noise_comps, var=out2$noise_var, ROI=out1$ROI_noise)
  )
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
#' @param ROI_noise List of numeric vectors, where each entry contains the label
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
#'  scores, one for each \code{ROI_noise}.
#'
#' @export
CompCor.HCP <- function(vol, labs=NULL, data_labs=NULL, data_cii=NULL, ROI_noise=NULL, noise_nPC=5, noise_erosion=0, verbose=TRUE){
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
    ROI_data <- array(labs %in% do.call(c, data_labs), dim=dim(labs))
    data <- list(matrix(vol[ROI_data], ncol=dim(vol)[4]))
  } else {
    if (is.character(data_cii)) { 
      if (requireNamespace("ciftiTools", quietly = TRUE)) {
        data_cii <- ciftiTools::read_cifti(data_cii) 
      } else {
        stop("Package `ciftiTools` required to read the CIFTI file. Please install it from the github repo `mandymejia/ciftiTools`.")
      }
    }
    data <- data_cii$data
    ROI_data <- NULL
  }

  # Check `ROI_noise`.
  if (verbose) { cat("Preparing noise ROIs (to use as regressors).\n") }
  ROI_noise.default <- list(
    wm_cort = c(3000:4035, 5001, 5002), 
    csf = c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255)
  )
  if (is.null(ROI_noise)) {
    ROI_noise <- ROI_noise.default
  } else if (is.character(ROI_noise)) {
    stopifnot(all(ROI_noise %in% c("wm_cort", "csf")))
    ROI_noise <- ROI_noise.default[unique(ROI_noise)]
  } else if (is.list(ROI_noise)) {
    ROI_noise <- ROI_noise
    stopifnot(is.numeric(do.call(c, data_labs)))
  } else {
    stop("`ROI_noise` argument is not in a recognized form.")
  }
  ROI_noise <- lapply(ROI_noise, function(x){array(labs %in% x, dim=dim(labs))})
  x <- CompCor(vol, ROI_noise=ROI_noise, noise_nPC=noise_nPC, noise_erosion=noise_erosion)

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

  list(data=data, ROI_data=ROI_data, noise=list(x$noise))
}