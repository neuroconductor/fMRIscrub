#' \code{clever} with CompCor for HCP data
#'
#' Wrapper to \code{\link{clever}} for HCP data & CompCor. The data from the 
#'  whole-brain NIFTI is used to obtain the noise ROIs, which are regressed
#'  from the greyordinate data in the CIFTI. 
#'
#' @param cii \code{"xifti"} (or file path to the CIFTI) from which the noise
#'  ROI components will be regressed. In the HCP, the corresponding file is e.g.
#'  "Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii"
#' @param nii \eqn{I \times J \times \K \times T} 
#'  NIFTI object or array (or file path to the NIFTI) which contains
#'  whole-brain data, including the noise ROIs. In the HCP, the corresponding
#'  file is e.g. "Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz"
#' @param nii_labels \eqn{I \times J \times K}
#'  NIFTI object or array (or file path to the NIFTI) which
#'  contains the corresponding labels to each voxel in \code{nii}. Values should
#'  be according to this table: 
#'  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT .
#'  In the HCP, the corresponding file is "ROIs/Atlas_wmparc.2.nii.gz". 
#' @param ROI_noise A list of numeric vectors. Each entry should represent labels
#'  in \code{nii_labels} belonging to a single noise ROI, named by that entry's 
#'  name. Or, this can be a character vector of at least one of the following: 
#'  \code{"wm_cort"} (cortical white matter), \code{"wm_cblm"} (cerebellar white
#'  matter), \code{"csf"} (cerebrospinal fluid). In the latter case, these labels
#'  will be used:
#'
#'  \describe{
#'    \item{wm_cort}{c(3000:4035, 5001, 5002)}
#'    \item{wm_cblm}{c(7, 46)}
#'    \item{csf}{c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255))}
#'  }
#'
#'  These default ROIs are based on this forum post: 
#'  https://www.mail-archive.com/hcp-users@humanconnectome.org/msg00931.html
#'
#'  Default: \code{c("wm_cort", "csf")}
#' @param noise_nPC,noise_erosion See \code{\link{clever}}.
#' @param ... Additional arguments to \code{\link{clever}}.
#'
#' @export
clever_HCP <- function(cii, nii, nii_labels, ROI_noise, noise_nPC=5, noise_erosion=NULL, ...){

  if (any(c("X", "ROI_data") %in% names(list(...)))) { 
    stop("`X` and `ROI_data` are overriden by the CIFTI/NIFTI file name\
    arguments for `clever_HCP`.") 
  }

  # `cii`
  if (is.character(cii)) { 
    if (requireNamespace("ciftiTools", quietly = TRUE)) {
      cii <- ciftiTools::read_cifti(cii) 
    } else {
      stop("Package `ciftiTools` required to read the CIFTI file. Please install it from the github repo `mandymejia/ciftiTools`.")
    }
  }
  stopifnot(all(names(cii) == c("data", "surf", "meta")))

  # `nii`
  if (is.character(nii)) {
    nii <- read_nifti(nii)
  }
  stopifnot(length(dim(nii))==4)

  # `labs`.
  if (is.character(nii_labels)) {
    nii_labels <- read_nifti(nii_labels)
  }
  stopifnot(length(dim(nii_labels))==3)
  stopifnot(dim(nii_labels)==dim(nii)[1:3])

  # Check `ROI_noise`.
  ROI_noise.default <- list(
    wm_cort = c(3000:4035, 5001, 5002), 
    wm_cblm =  c(7, 46),
    csf = c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255)
  )
  if (is.null(ROI_noise)) {
    ROI_noise <- ROI_noise.default[c("wm_cort", "csf")]
  } else if (is.character(ROI_noise)) {
    stopifnot(all(ROI_noise %in% names(ROI_noise.default)))
    ROI_noise <- ROI_noise.default[unique(ROI_noise)]
  } else if (is.list(ROI_noise)) {
    ROI_noise <- ROI_noise
    stopifnot(is.numeric(do.call(c, ROI_noise)))
  } else {
    stop("`ROI_noise` argument is not in a recognized form.")
  }
  ROI_noise <- lapply(ROI_noise, function(x){array(nii_labels %in% x, dim=dim(nii_labels))})

  temp <- format_data(
    nii, ROI_data=NULL, ROI_noise = ROI_noise, 
    noise_erosion=noise_erosion, noise_nPC=noise_nPC
  )

  out <- clever(
    t(do.call(rbind, cii$data)), 
    ROI_noise=temp$X_noise, noise_nPC=as.numeric(temp$noise_nPC), 
    ...
  )
  out$ROIs <- c(out$ROIs, temp$ROI_noise)
  out
}