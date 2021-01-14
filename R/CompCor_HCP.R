#' Run CompCor on the HCP data
#'
#' Wrapper to \code{\link{clever_multi}} for computing CompCor (and other outlyingness
#'  measures) on HCP data. The whole-brain NIFTI is used to obtain the noise
#'  ROIs, which are regressed from the greyordinate data in the CIFTI. 
#'
#' @param nii \eqn{I \times J \times \K \times T} 
#'  NIFTI object or array (or file path to the NIFTI) which contains
#'  whole-brain data, including the noise ROIs. In the HCP, the corresponding
#'  file is e.g. "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz"
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
#' @param noise_erosion  The number of voxel layers to erode the noise ROIs by. 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its 
#'  elements will be recycled. Default: \code{NULL}, which will use a value of
#'  0 (do not erode the noise ROIs).
#' @param brainstructures Choose among "left", "right", and "subcortical".
#'  Default: \code{c("left", "right")} (cortical data only)
#' @param cii \code{"xifti"} (or file path to the CIFTI) from which the noise
#'  ROI components will be regressed. In the HCP, the corresponding file is e.g.
#'  "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii"
#' @param frames A numeric vector indicating the timepoints to compute 
#'  CompCor for, or \code{NULL} (default) to use all frames. (Indexing begins
#'  with 1, so the first timepoint has index 1 and the last has the same index
#'  as the length of the scan.)
#' @param center,scale Center the columns of the data by median, and scale the
#'  columns of the data by MAD? Default: \code{TRUE} for both. Affects both
#'  \code{X} and the noise data.
#' @param DCT Detrend the columns of the data using the discrete cosine
#'  transform (DCT)? Use an integer to indicate the number of cosine bases to 
#'  use for detrending. Use \code{0} (default) to forgo detrending. 
#' 
#'  The data must be centered, either before input or with \code{center}.
#' @param nuisance_too A matrix of nuisance signals to regress from the data
#'  before, i.e. a "design matrix." Should have \eqn{T} rows. Nuisance
#'  regression will be performed simultaneously with DCT detrending if 
#'  applicable. \code{NULL} to not add additional nuisance regressors. Affects 
#'  both \code{X} and the noise data.
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @export
CompCor_HCP <- function(
  nii, nii_labels, 
  ROI_noise=c("wm_cort", "csf"), noise_nPC=5, noise_erosion=NULL, 
  frames=NULL, cii=NULL, brainstructures=c("left", "right"),
  center = TRUE, scale = TRUE, DCT = 0, nuisance_too = NULL,
  verbose=FALSE){

  # `nii`
  if (is.character(nii)) {
    cat("Reading data NIFTI.\n")
    nii <- read_nifti(nii)
  }
  stopifnot(length(dim(nii))==4)

  # `labs`.
  if (is.character(nii_labels)) {
    cat("Reading labels NIFTI.\n")
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
  ROI_noise <- lapply(
    ROI_noise, 
    function(x){array(nii_labels %in% x, dim=dim(nii_labels))}
  )

  T_ <- dim(nii)[4]
  if (!is.null(frames)) {
    stopifnot(all(frames %in% seq(T_)))
    nii <- nii[,,,frames]
    if (length(frames) < 10) {
      warning("There are very few frames.\n")
    }
  }

  cat("Computing noise components.\n")
  out <- CompCor(
    nii, ROI_data=NULL, ROI_noise=ROI_noise, 
    noise_erosion=noise_erosion, noise_nPC=noise_nPC,
    center=center, scale=scale, DCT=DCT, nuisance_too=nuisance_too
  )$noise

  # `cii`
  if (!is.null(cii)) {
    if (is.character(cii)) { 
      if (requireNamespace("ciftiTools", quietly = TRUE)) {
        cii <- ciftiTools::read_cifti(cii, brainstructures=brainstructures, verbose=verbose) 
      } else {
        stop("Package `ciftiTools` required to read the CIFTI file. Please install\
        it from the github repo `mandymejia/ciftiTools`.")
      }
    }
    stopifnot(all(names(cii) == c("data", "surf", "meta")))

    cii <- do.call(rbind, cii$data)

    if (!is.null(frames)) {
      stopifnot(all(frames %in% seq(ncol(cii))))
      cii <- cii[,frames]
    }

    out$data <- CompCor.regress(cii, out$PCs)
  }

  out
}