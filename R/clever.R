#' Identify outliers with \code{clever}
#' 
#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' \code{clever} will use all combinations of the requested projections and 
#'  outlyingness measures which make sense. For example, if  
#'  \code{projection=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{out_meas=c("leverage", "robdist")} then these five
#'  combinations will be used: PCATF with leverage, PCA + variance with 
#'  leverage, PCA + variance with robust distance, PCA + kurtosis with leverage,
#'  and PCA + kurtosis with robust distance. Each method combination will yield 
#'  its own results.
#' 
#' Each voxel timecourse is centered on its median and scaled by 1.4826 times
#'  the median of the absolute values (a robust measure of standard deviation).
#'  This differs from the centering method of DVARS (Afyouni and Nichols, 2018),
#'  so results may differ.
#' 
#' @param X Wide numerical data matrix (\eqn{T observations \times V variables}, \eqn{T << V}).
#'  For example, if \code{X} represents an fMRI run, \eqn{T} should be the number
#'  of timepoints and \eqn{V} should be the number of brainordinate vertices/voxels.
#' 
#'  Or, a 4D array or NIFTI or file path to a NIFTI (\eqn{I \times J \times K \times T} 
#'  observations), in which case \code{ROI_data} must be provided.
#' 
#'  Or, a \code{ciftiTools} code{"xifti"} object or a file path to a CIFTI
#'  (\eqn{V_{left} + V_{right} + V_{subcortical} \times T observations}).
#' @param measures Character vector indicating the outlyingness measures to 
#'  compute. Choose at least one of the following: 
#' 
#'  \describe{
#'    \item{\code{"leverage"}}{PCA leverage (the mean of the squared PC scores)}
#'    \item{\code{"robdist"}}{Robust Mahalanobis-based distance}
#'    \item{\code{"CompCor"}}{anatomical CompCor based on the ROIs. \code{X} must
#'      be a 4D array or NIFTI, and \code{ROI_data} and \code{ROI_noise} are required.}
#'    \item{\code{"DVARS"}}{traditional-DVARS, as well as Delta-percent-DVARS and
#'      z-score-DVARS (Afyouni and Nichols, 2018)}
#'    \item{\code{"FD"}}{Framewise Displacement. \code{X_motion} is required.}
#'  }
#' 
#'  Use \code{"all"} to select all outlyingness measures. Default: \code{"leverage", "DVARS"}.
#' @param ROI_data Indicates the data ROI. 
#' 
#'  If \code{X} is a matrix, this must be a length \eqn{V} logical vector, where
#'  the data ROI is indicated by \code{TRUE} values. If not provided, all 
#'  columns of \code{X} will be included in the data ROI (all \code{TRUE}).
#' 
#'  If \code{X} is an array or NIFTI, this must be either a vector of values
#'  to expect for out-of-mask voxels in \code{X}, or a (file path to a) 3D NIFTI.
#'  In the latter case, each of the volume dimensions should match the first
#'  three dimension of \code{X}. Voxels in the data ROI should be indicated by
#'  \code{TRUE} and all other voxels by \code{FALSE}. If not provided,
#'  will be set to \code{c(0, NA, NaN)} (include all voxels which are not constant
#'  \code{0}, \code{NA}, or \code{NaN}).
#' 
#'  If \code{X} is a \code{"xifti"} this should not be used. All data locations
#'  will be used. 
#' @param ROI_noise Indicates the noise ROI. Only used if the \code{"CompCor"} 
#'  measure is requested.
#'  
#'  If \code{X} is a matrix, this must be a list of length \eqn{V} logical
#'  vectors, or a list of matrices with \code{T} rows. The names of each entry should
#'  indicate the name of the noise ROI, e.g. \code{"white_matter"} and \code{"csf"}.
#'  In the first case, \code{TRUE} values should indicate the locations of \code{X} 
#'  within that noise ROI. Since the ROIs must not overlap, the masks must be 
#'  mutually exclusive with each other, and with \code{ROI_data}. In the second
#'  case, the rows of the matrix must represent noise brainordinate timecourses,
#'  separate from \code{X}. 
#' 
#'  If \code{X} is an array or NIFTI, this must be a list of (file paths to) 3D 
#'  NIFTIs or arrays, or a list of matrices with \code{T} rows. The names of 
#'  each entry should indicate the name of the noise ROI, e.g. 
#'  \code{"white_matter"} and \code{"csf"}. In the first case, each of the volume 
#'  dimensions should match the first three dimensions of \code{X}. Voxels in 
#'  each noise ROI should be indicated by \code{TRUE} and all other voxels by 
#'  \code{FALSE}. Since the ROIs must not overlap, the masks must be mutually 
#'  exclusive with each other, and with \code{ROI_data}. In the second case,
#'  the rows of the matrix must represent noise brainordinate timecourses,
#'  separate from \code{X}. 
#' 
#'  If \code{X} is a \code{"xifti"}, this must be a list of matrices with 
#'  \code{T} rows. The names of each entry should indicate the name of the noise
#'  ROI, e.g. \code{"white_matter"} and \code{"csf"}. The rows of the matrix 
#'  must represent noise brainordinate timecourses, separate from \code{X}. 
#' @param X_motion Only used if the \code{"FD"} measure is requested. An 
#'  \eqn{N \times 6} matrix in which the first three columns represent the
#'  translational realignment parameters (mm), and the second three columns represent
#'  the rotational realignment parameters in (radians). To convert radians to mm,
#'  the displacement on a sphere of radius 50 mm will be computed.
#' @param projections Only applies to the \code{"leverage"} and \code{"robdist"} 
#'  measures. These work by projecting the data onto directions likely to 
#'  contain outlier information. Choose at least one of the following:
#' 
#'  \describe{
#'    \item{\code{"PCA_var"}}{PCA using the PCs of above-average variance. Compatible with both leverage and robust distance.}
#'    \item{\code{"PCA_kurt"}}{PCA using the PCs of high kurtosis and above-average variance. Compatible with both leverage and robust distance.}
#'    \item{\code{"PCATF"}}{PCATF using the trend-filtered PCs of above-average variance. Compatible with only leverage.}
#'  }
#'  
#'  Use \code{"all"} to use all projection methods. Default: \code{"PCA_kurt"}.
#'  
#'  Each compatible measure + projection combination will yield its own result.
#'  For example, if all projections are used and both leverage and robust distance
#'  are requested, the results will be: leverage of high-variance PCs, leverage
#'  of high-kurtosis PCs, leverage of trend-filtered PCs, robust distance of
#'  high-variance PCs, and robust distance of high-kurtosis PCs.
#' @param compute_PC_dirs Only applies to the \code{"leverage"} and \code{"robdist"} 
#'  measures. Should the PCA and PCATF principal directions be computed? 
#'  Default: \code{FALSE} (conserves memory). Required to use \code{\link{leverage_images}}.
#' @param detrend Only applies to the \code{"leverage"} and \code{"robdist"} 
#'  measures. Detrend the PCs before measuring kurtosis and before measuring
#'  leverage or robust distance? Default: \code{TRUE}.
#' 
#'  Detrending is highly recommended for time-series data, especially if there 
#'  are many time points or evolving circumstances affecting the data. Additionally,
#'  for the kurtosis-based projection, trends can induce positive or negative kurtosis,
#'  contaminating the connection between high kurtosis and outlier presence. 
#'  
#'  Detrending should not be used with non-time-series data because the 
#'  observations are not temporally related.
#' 
#'  In addition to \code{TRUE} and \code{FALSE}, a third option \code{"kurtosis"}
#'  can be used to only detrend the PCs for the purpose of measuring kurtosis, 
#'  and not for the actual outlyingness measurement.
#' 
#'  This option will not affect the PCATF PCs, which are never detrended.
#' @param noise_erosion Only applies to the \code{"CompCor"} measure. The number
#'  of voxel layers to erode the noise ROIs by. Should be a list or numeric
#'  vector with the same length as \code{ROI_noise}. It will be matched to each
#'  ROI based on the name of each entry, or if the names do not match or are missing,
#'  the order of entries. If it is an unnamed vector, its elements will be recycled.
#'  Default: \code{NULL}, which will use a value of 0 (do not erode the noise ROIs).
#' @param noise_nPC Only applies to the \code{"CompCor"} measure. The number of 
#'  principal components to compute for each noise ROI. Should be a list or numeric
#'  vector with the same length as \code{ROI_noise}. It will be matched to each
#'  ROI based on the name of each entry, or if the names do not match or are missing,
#'  the order of entries. If it is an unnamed vector, its elements will be recycled.
#'  Default: \code{5} (compute the top 5 PCs for each noise ROI).
#' @param PCATF_kwargs Options for the \code{"PCATF"} projection. Valid entries 
#'  are: 
#'  
#'  \describe{
#'    \item{K}{Maximum number of PCs to compute. Default: \code{100}. Cannot be set above 100.}
#'    \item{lambda}{Trend-filtering parameter. Default: \code{0.05}.}
#'    \item{niter_max}{Maximum number of iterations. Default: \code{1000}.}
#'    \item{verbose}{Print updates? Default: \code{FALSE}.}
#'  }
#' @param kurt_quantile Only applies to \code{"PCA_kurt"} projection. 
#'  What cutoff quantile for kurtosis should be used to select the PCs? 
#'  Default: \code{0.95}.
#' @param get_outliers Should outliers be flagged based on cutoffs? Default: \code{TRUE}.
#' @param outlier_cutoffs Named list indicating the cutoff for each outlyingness measure.
#'  Only used if \code{get_outliers}. Each cutoff is specified in a particular way:
#' 
#'  \describe{
#'    \item{\code{"leverage"}}{Minimum leverage value, in multiples of the median leverage. Default: \code{4} (will flag leverage scores more than four times the median).}
#'    \item{\code{"robdist"}}{Minimum robust distance quantile, based on the estimated F distribution. Default: \code{.9999} (will flag robust distance scores above the 99.99th percentile).}
#'    \item{\code{"CompCor"}}{Not applicable: not used for outlier detection. Should not be specified.}
#'    \item{\code{"DVARS"}}{Minimum traditional-DVARS value, or "Afyouni_Nichols" (default), in which case
#'      their dual cutoff will be used (Delta-percent-DVARS of more than +5\% or 
#'      less than -5\% AND z-score=DVARS greater than the right-tail
#'      5\% significance level with Bonferroni FWER correction).}
#'    \item{\code{"FD"}}{Minimum FD. Default: 0.5 mm (will flag FD scores greater than 0.5 mm).}
#'  }
# @param R_true The N x N correlation matrix, if known. Used for the bootstrap
#  robust distance method.
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @return A \code{"clever"} object, i.e. a list with components
#' \describe{
#'  \item{measures}{
#'    A list of the outlyingness measures (only those requested will be included):
#'    \describe{
#'      \item{leverage__PCA_var}{Leverage values of the "PCA_var" projection.}
#'      \item{leverage__PCA_kurt}{Leverage values of the "PCA_kurt" projection.}
#'      \item{leverage__PCATF}{Leverage values of the "PCATF" projection.}
#'      \item{robdist__PCA_var}{Robust distance values of the "PCA_var" projection.}
#'      \item{robdist__PCA_kurt}{Robust distance values of the "PCA_kurt" projection.}
#'      \item{CompCor_[Noise1]_PC1}{First PC of the first noise ROI.}
#'      \item{...}{...}
#'      \item{CompCor_[Noisek]_PCn}{nth PC of the kth (last) noise ROI.}
#'      \item{DVARS}{Traditional DVARS values.}
#'      \item{DVARS_DPD}{Delta-percent-DVARS values.}
#'      \item{DVARS_ZD}{z-score-DVARS values.}
#'      \item{FD}{Framewise Displacement values}
#'    }
#'  }
#'  \item{outlier_cutoffs}{
#'    A list of the outlier cutoffs for each measure (see the \code{outlier_cutoffs}
#'    argument; only those requested will be included):
#'    \describe{
#'      \item{leverage__PCA_var}{Minimum leverage.}
#'      \item{leverage__PCA_kurt}{Minimum leverage.}
#'      \item{leverage__PCATF}{Minimum leverage.}
#'      \item{robdist__PCA_var}{Minimum robust distance.}
#'      \item{robdist__PCA_kurt}{Minimum robust distance.}
#'      \item{DVARS}{Minimum DVARS.}
#'      \item{DVARS_DPD}{Minimum absolute value of the Delta-percent-DVARS.}
#'      \item{DVARS_ZD}{Minimum z-score-DVARS.}
#'      \item{FD}{Minimum Framewise Displacement.}
#'    }
#'  }
#'  \item{outlier_flags}{
#'    Applies \code{outlier_cutoffs} to \code{measures}: each is a logical vector
#'    with \code{TRUE} values indicating suspected outlier presence. 
#'  }
#'  \item{ROIs}{
#'    \describe{
#'      \item{data}{The mask of locations in the data ROI.}
#'      \item{[Noise1]}{The mask of locations in the first noise ROI, if it was relative to \code{X}.}
#'      \item{...}{...}
#'      \item{[Noisek]}{The mask of locations in the kth (last) noise ROI, if it was relative to \code{X}.}
#'    }
#'  }
#'  \item{PCA}{
#'    If the "PCA_var" or "PCA_kurt" projections were used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{N x Q} PC score matrix. Only PCs with above-average variance will be included.}
#'      \item{U_dt}{The \eqn{N x Q} detrended PC score matrix. Included only if \code{detrend}}
#'      \item{D}{The variance of each PC. Only PCs with above-average variance will be included.}
#'      \item{V}{The \eqn{P x Q} PC directions matrix. Included only if \code{!compute_PC_dirs}}
#'      \item{kurt_idx}{The length \code{Q} kurtosis rankings, with 1 indicating the highest-kurtosis PC 
#'        (among those of above-average variance) and \code{NA} indicating a PC with kurtosis below
#'        the quantile cutoff. Only included if the "PCA_kurt" projection was used.}
#'    }
#'  }
#'  \item{PCATF}{
#'    If the "PCATF" projection was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{N x Q} PC score matrix. Only PCs with above-average variance will be included.}
#'      \item{D}{The variance of each PC. Only PCs with above-average variance will be included.}
#'      \item{V}{The \eqn{P x Q} PC directions matrix. Included only if \code{!compute_PC_dirs}}
#'    }
#'  }
#'  \item{robdist}{
#'    If the "robdist" method was used, this will be a list with components:
#'    \describe{
#'      \item{PCA_var}{
#'      If the "PCA_var" projection was used, this will be a list with components:
#'        \describe{
#'          \item{inMCD}{Logical vector indicating whether each observation was in the MCD estimate.}
#'          \item{outMCD_scale}{The scale for out-of-MCD observations.}
#'          \item{Fparam}{Named numeric vector: \code{c}, \code{m}, \code{df1}, and \code{df2}.}
#'        }
#'      }
#'      \item{PCA_kurt}{
#'      If the "PCA_kurt" projection was used, this will be a list with components:
#'        \describe{
#'          \item{inMCD}{Logical vector indicating whether each observation was in the MCD estimate.}
#'          \item{outMCD_scale}{The scale for out-of-MCD observations.}
#'          \item{Fparam}{Named numeric vector: \code{c}, \code{m}, \code{df1}, and \code{df2}.}
#'        }
#'      }
#'    }
#'  }
#' }
#'
#' @importFrom robustbase rowMedians
#' @import stats
#'
#' @export
#'
#' @examples
#' n_voxels = 1e4
#' n_timepoints = 100
#' X = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' clev = clever(X)
clever = function(
  X,
  measures=c("leverage", "DVARS"),
  ROI_data=c(NA, NaN), ROI_noise=NULL, X_motion=NULL,
  projections = "PCA_kurt", compute_PC_dirs=FALSE,
  detrend=TRUE,
  noise_nPC=5, noise_erosion=NULL,
  PCATF_kwargs=NULL, kurt_quantile=.95,
  get_outliers=TRUE, 
  outlier_cutoffs=list(leverage=4, robdist=.9999, DVARS="Afyouni_Nichols", FD=.5),
  verbose=FALSE){

  # Define the cutoff value for detecting zero variance/MAD voxels
  TOL <- 1e-8

  # ----------------------------------------------------------------------------
  # Check arguments. -----------------------------------------------------------
  # ----------------------------------------------------------------------------

  valid_measures <- c("leverage", "robdist", "CompCor", "DVARS", "FD") # robdist_bootstrap
  valid_projections <- c("PCA_var", "PCA_kurt", "PCATF")

  if ("all" %in% measures) {
    measures <- valid_measures
  } else {
    measures <- unique(match.arg(measures, valid_measures, several.ok=TRUE))
  }

  if ("all" %in% projections) {
    projections <- valid_projections
  } else {
    projections <- unique(match.arg(projections, valid_projections, several.ok=TRUE))
  }

  use_PCA <- FALSE
  use_PCATF <- FALSE
  if ("leverage" %in% measures) {
    if ("PCA_var" %in% projections | "PCA_kurt" %in% projections) { use_PCA <- TRUE }
    if ("PCATF" %in% projections) { use_PCATF <- TRUE }
    measures <- measures[measures != "leverage"]
    measures <- c(measures, paste0(leverage, "__", projections))
  }
  if ("robdist" %in% measures) {
    if ("PCA_var" %in% projections | "PCA_kurt" %in% projections) { use_PCA <- TRUE }
    measures <- measures[measures != "robdist"]
    measures <- c(measures, paste0(robdist, "__", projections))
    measures <- measures[measures != "robdist__PCATF"] # not compatible
  }

  if(identical(projection, "all")){
    projection <- all_projection
  } else {
    projection <- match.arg(projection, all_projection, several.ok=TRUE)
  }

  if(identical(out_meas, "all")){
    out_meas <- all_out_meas
  } else {
    out_meas <- match.arg(out_meas, all_out_meas, several.ok=TRUE)
  }

  is.TRUEorFALSE <- function(x) { length(x)==1 && is.logical(x) }
  stopifnot(is.TRUEorFALSE(compute_PC_dirs))
  stopifnot(is.TRUEorFALSE(get_outliers))
  stopifnot(is.TRUEorFALSE(verbose))

  stopifnot(length(detrend)==1)
  detrend <- switch(as.character(detrend),
    `TRUE` = c("PCs", "kurtosis"),
    kurtosis = "kurtosis",
    `FALSE` = NULL
  )

  if(!identical(PCATF_kwargs, NULL)){
    names(PCATF_kwargs) <- match.arg(
      names(PCATF_kwargs), c("K", "lambda", "niter_max", "TOL", "verbose"),
      several.ok=TRUE)
    if(length(PCATF_kwargs) != length(unique(unlist(PCATF_kwargs)))){
      stop("Duplicate PCATF_kwargs were given.")
    }
  }

  stopifnot((kurt_quantile < 1) & (kurt_quantile > 0))

  outlier_cutoffs <- as.list(outlier_cutoffs)
  outlier_cutoffs_defaults <- list(leverage=4, robdist=.9999, DVARS="Afyouni_Nichols", FD=.5)
  for (meas in names(outlier_cutoffs_defaults)) {
    if (!(meas %in% names(outlier_cutoffs_defaults))) { 
      outlier_cutoffs[[meas]] <- outlier_cutoffs_defaults[[meas]]
    }
  }
  stopifnot(all(sapply(outlier_cutoffs, length) == 1))
  stopifnot(is.numeric(outlier_cutoffs$leverage) & outlier_cutoffs$leverage > 0)
  stopifnot(is.numeric(outlier_cutoffs$robdist))
  stopifnot((outlier_cutoffs$robdist > 0) & (outlier_cutoffs$robdist < 1))
  stopifnot(is.numeric(outlier_cutoffs$DVARS) | outlier_cutoffs_defaults$DVARS=="Afyouni_Nichols")
  stopifnot(is.numeric(outlier_cutoffs$FD) & outlier_cutoffs$FD > 0)

  # ----------------------------------------------------------------------------
  # Check data. ----------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # X
  if (is.matrix(X)) {
    Vpre_ <- ncol(X); T_ <- nrow(X)
    if(Vpre_ < T_){
      warning(
        "Data matrix has more rows than columns. Check that observations\
        are in rows and variables are in columns."
      )
    }
    X_type <- "vectorized"
  } else if (is.array(X) && length(dim(X))==4) {
    T_ <- dim(X)[4]
    X_type <- "volume"
  } else if (is.character(X)) {
    if (endsWith(X, ".dtseries.nii") | endsWith(X, ".dscalar.nii")) {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read `X`. Please install it", call. = FALSE)
      }
      X_cifti <- read_cifti(X, brainstructures="all")
      X <- t(do.call(rbind, X_cifti$data))
      Vpre_ <- ncol(X); T_ <- nrow(X)
      X_type <- "cifti"
    } else {
      X <- read_nifti(X)
      X_type <- "volume"
    }
  } else if (inherits(X, "xifti")) {
    X_cifti <- X
    X <- t(do.call(rbind, X_cifti$data))
    Vpre_ <- ncol(X); T_ <- nrow(X)
    X_type <- "cifti"
  } else {
    stop("`X` must be a matrix, array, NIFTI, path to a NIFTI, CIFTI, or path to a CIFTI.")
  }

  # ROI_data
  ROI_data_was_null <- is.null(ROI_data)
  if (X_type == "vectorized") {
    if (ROI_data_was_null) { ROI_data <- rep(TRUE, Vpre_) }
    ROI_data <- as.vector(ROI_data)
    stopifnot(length(ROI_data) == Vpre_)
    ROI_data <- as.logical(ROI_data)
  } else if (X_type == "volume") {
    if (ROI_data_was_null) { ROI_data <- c(0, NA, NaN) }
    if (is.character(ROI_data)) {
      ROI_data <- read_nifti(ROI_data)
    }
    if (is.vector(ROI_data)) {
      ROI_data <- X
      ROI_data[,,] <- X %in% ROI_data
    } else if (is.array(ROI_data)) {
      stopifnot(all(dim(ROI_data) == dim(X)[1:3]))
      ROI_data <- ROI_data
      ROI_data[,,] <- as.logical(ROI_data)
    }
  } else if (X_type == "xifti") {
    ROI_data <- rep(TRUE, Vpre_)
  } else { stop("Internal error: unrecognized `X_type`") }

  # ROI_noise
  if ("CompCor" %in% measures) {
    if (is.null(ROI_noise)) { stop("`CompCor` requires the noise ROIs.") }
    if (!is.list(ROI_noise)) { ROI_noise <- list(Noise1=ROI_noise) }
    if (is.null(names(ROI_noise))) { names(ROI_noise) <- paste0("Noise", 1:length(ROI_noise)) }
    stopifnot(length(names(ROI_noise)) == length(unique(names(ROI_noise))))

    # noise_nPC
    noise_nPC <- as.list(noise_nPC)
    if (is.null(names(noise_nPC))) {
      noise_nPC <- noise_nPC[rep(1:length(noise_nPC), length(ROI_noise))[1:length(ROI_noise)]]
      names(noise_nPC) <- names(ROI_noise)
    } else {
      stopifnot(all(sorted(names(noise_nPC)) == sorted(names(ROI_noise))))
    }  

    # noise_erosion
    if (is.null(noise_erosion)) { 
      noise_erosion = 0
    } else {
      if (!any(sapply(ROI_noise, is.array))) {
        warning("`noise_erosion` was provided, but there are no array/NIFTI noise ROIs to erode.")
      }
    }
    noise_erosion <- as.list(noise_erosion)
    if (is.null(names(noise_erosion))) {
      noise_erosion <- noise_erosion[rep(1:length(noise_erosion), length(ROI_noise))[1:length(ROI_noise)]]
      names(noise_erosion) <- names(ROI_noise)
    } else {
      stopifnot(all(sorted(names(noise_erosion)) == sorted(names(ROI_noise))))
    }

    X_noise <- vector("list", length(ROI_noise)); names(X_noise) <- names(ROI_noise)
    for (ii in 1:length(ROI_noise)) {
      if (is.null(ROI_noise[[ii]])) { ROI_noise[[ii]] <- NULL; next }
      if (X_type == "vectorized") {
        if (is.vector(ROI_noise[[ii]])) {
          stopifnot(length(ROI_noise[[ii]]) == Vpre_)
          ROI_noise[[ii]] <- as.logical(ROI_noise[[ii]])
          X_noise[[ii]] <- X[,ROI_noise[[ii]]]
        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
        } else {
          stop("Each entry in `ROI_noise` must be a logical vector, or matrix with the same number of rows as `X`.")
        }
      } else if (X_type == "volume") {
        if (is.character(ROI_noise[[ii]])) {
          if (!file.exists(ROI_noise[[ii]])) { stop(paste("The `ROI_noise` entry", ROI_noise[[ii]], "is not an existing file.")) }
          ROI_noise[[ii]] <- read_nifti(ROI_noise[[ii]])
        }
        if (is.array(ROI_noise[[ii]])) {
          stopifnot(all(dim(ROI_noise[[ii]]) == dim(X)[1:3]))
          ROI_noise[[ii]][,,] <- as.logical(ROI_noise[[ii]])
          X_noise[[ii]] <- t(matrix(X[ROI_noise[[ii]]],ncol=dim(X)[4]))
        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
        } else {
          stop("Each entry in `ROI_noise` must be a logical array, or matrix with the same number of rows as `X`.") 
        }
      } else if (X_type == "xifti") {
        stopifnot(is.matrix(ROI_noise[[ii]]))
        stopifnot(nrow(ROI_noise[[ii]]) == T_)
        X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[[ii]] <- NULL
      } else { stop("Internal error: unrecognized `X_type`") }
    }
  }

  if (!all(sapply(ROI_noise, is.null))) {
    # check that ROI are mutually exclusive
    all_noise_ROIs <- apply(do.call(rbind, ROI_noise), 2, sum)
    if (!all(all_noise_ROIs < 2)) {
      stop("The noise ROIs must all be mutually exclusive.")
    }
    all_noise_ROIs <- all_noise_ROIs > 0
    if (ROI_data_was_null) { 
      if (X_type == "matrix") {
        ROI_data[all_noise_ROIs] <- FALSE
      } else if (X_type == "array") {
        ROI_data[,,][all_noise_ROIs] <- FALSE
      }
    } else {
      if (any(all_noise_ROIs & as.vector(ROI_data))) {
        stop("The noise ROIs must not overlap with the data ROI.")
      }
    }
  }

  if (X_type == "array") {
    X <- t(matrix(X[!all_noise_ROIs], ncol=dim(X)[4]))
  }

  # make sure nPCs greater than the rank of the noise ROIs!
  # similarly for data
  # ...

  # Get data dimensions.
  Vpre_ <- ncol(X); T_ <- nrow(X)
  if(Vpre_ < T_){
    warning(
      "Data matrix has more rows than columns. Check that observations\
      are in rows and variables are in columns."
    )
  }

  return(list(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise
  ))

  # Collect all the methods to compute.
  methods <- all_valid_methods[
    all_valid_methods %in% outer(projection, out_meas, paste, sep='__')
  ]
  if(length(methods) < 1){
    stop(
      "No valid method combinations. Check that `projection` and\
      `out_meas` are compatible.\n"
    )
  }
  outlier_measures <- outlier_lev_imgs <- setNames(vector("list", length(methods)), methods)
  if(get_outliers){
    outlier_cutoffs <- outlier_flags <- setNames(vector("list", length(methods)), methods)
  }
  robdist_info <- vector("list")

  # ----------------------------------------------------------------------------
  # Center and scale the data. -------------------------------------------------
  # Do it here instead of calling `scale_med` to save memory. ------------------
  # ----------------------------------------------------------------------------

  if(verbose){ cat("Centering and scaling the data matrix.\n") }
  # Transpose.
  X <- t(X)
  #	Center.
  X <- X - c(rowMedians(X, na.rm=TRUE))
  # Scale.
  mad <- 1.4826 * rowMedians(abs(X), na.rm=TRUE)
  const_mask <- mad < TOL
  if(any(const_mask)){
    if(all(const_mask)){
    stop("All voxels are zero-variance.\n")
    } else {
      warning(paste0("Warning: ", sum(const_mask),
      " constant voxels (out of ", length(const_mask),
      "). These will be removed for estimation of the covariance.\n"))
    }
  }
  mad <- mad[!const_mask]
  X <- X[!const_mask,]
  X <- X/c(mad)
  # Revert transpose.
  X <- t(X)
  N_ <- ncol(X)

  # ----------------------------------------------------------------------------
  # Compute DVARS. -------------------------------------------------------------
  # ----------------------------------------------------------------------------

  if(DVARS){
    if(verbose){ cat("Computing DVARS.\n") }
    X_DVARS <- DVARS(X, normalize=FALSE, norm_I=100, verbose=verbose)
    
    outlier_measures$DVARS_DPD <- X_DVARS$DPD
    outlier_measures$DVARS_ZD <- X_DVARS$ZD

    if(get_outliers){
      outlier_cutoffs$DVARS_DPD <- 5
      outlier_cutoffs$DVARS_ZD <- qnorm(1-.05/T_)
      
      outlier_flags$DVARS_DPD <- X_DVARS$DPD > outlier_cutoffs$DVARS_DPD
      outlier_flags$DVARS_ZD <- X_DVARS$ZD > outlier_cutoffs$DVARS_ZD
      outlier_flags$DVARS_both <- outlier_flags$DVARS_DPD & outlier_flags$DVARS_ZD
    }
  }

  # ----------------------------------------------------------------------------
  # Make projections. ----------------------------------------------------------
  # ----------------------------------------------------------------------------

  # Compute the PC scores (and directions, if leverage images or PCATF are desired).
  solve_dirs <- (lev_images) || ("PCATF" %in% projection)
  if(verbose){
    cat(paste0(
      "Computing the",
      ifelse(
        ("PCA_var" %in% projection) | ("PCA_kurt" %in% projection),
        ifelse("PCATF" %in% projection, " normal and trend-filtered", ""),
        ifelse("PCATF" %in% projection, " trend-filtered", "INTERNAL ERROR")
      ),
      " PC scores",
      ifelse(solve_dirs, " and directions", ""), ".\n"
    ))
  }
  if(solve_dirs){
    X.svd <- svd(X)
    if(!("PCATF" %in% projection)){ rm(X) }
  } else {
    # Conserve memory by using `XXt`
    XXt <- tcrossprod(X)
    if(!("PCATF" %in% projection)){ rm(X) }
    X.svd <- svd(XXt)
    rm(XXt)
    X.svd$d <- sqrt(X.svd$d)
    X.svd$v <- NULL
  }

  # Compute PCATF, if requested.
  if("PCATF" %in% projection){
    X.svdtf <- do.call(
      PCATF, 
      c(list(X=X, X.svd=X.svd, solve_directions=solve_dirs), PCATF_kwargs)
    )
    # The PC directions were needed to compute PCATF. If leverage images 
    #   are not wanted, we can now delete the directions to save space.
    if(!lev_images){ X.svd$v <- NULL }

    # Remove trend-filtered PCs with constant scores.
    # TO DO: Reconsider adding this back?
    tf_zero_var <- apply(X.svdtf$u, 2, var) < TOL
    if(any(tf_zero_var)){
      if(all(tf_zero_var)){
        stop("Error: All trend-filtered PC scores are zero-variance.")
      }
      # warning(paste(
      #   "Warning:", sum(tf_zero_var), 
      #   "trend-filtered PC scores are zero-variance. Removing these PCs.\n"
      # ))
      # X.svdtf$u <- X.svdtf$u[,!tf_zero_var]
      # X.svdtf$d <- X.svdtf$d[!tf_zero_var]
      # if(lev_images){ X.svdtf$v <- X.svdtf$v[,!tf_zero_var] }
    }
  }
  gc()

  # Choose which PCs to retain for each projection.
  projection <- setNames(vector("list", length(projection)), projection)
  for (ii in 1:length(projection)) {
    proj_ii_name <- names(projection)[ii]

    if(verbose){
      cat(switch(proj_ii_name,
        PCA_var = "Identifying the PCs with high varaince.\n",
        PCA_kurt = "Identifying the PCs with high kurtosis.\n",
        PCATF = "Identifying the trend-filtered PCs with high varaince.\n"
      ))
    }

    # Identify the indices to keep.
    chosen_PCs <- switch(proj_ii_name,
      PCATF = 1:ncol(X.svdtf$u),
      PCA_var = choose_PCs.variance(svd=X.svd),
      PCA_kurt = choose_PCs.kurtosis(
        svd=X.svd, 
        kurt_quantile=kurt_quantile, detrend="kurtosis" %in% detrend
      )
    )
    ## kurtosis order =/= index order
    chosen_PCs_ordered <- chosen_PCs[order(chosen_PCs)]
    
    # Get the PC subset. 
    proj_ii_svd <- switch(proj_ii_name,
      PCATF = X.svdtf,
      PCA_var = X.svd,
      PCA_kurt = X.svd
    )
    proj_ii_svd$u <- proj_ii_svd$u[,chosen_PCs_ordered,drop=FALSE]
    proj_ii_svd$d <- proj_ii_svd$d[chosen_PCs_ordered]
    if("PCs" %in% detrend && proj_ii_name != "PCATF"){
      proj_ii_svd$u_detrended <- proj_ii_svd$u - apply(proj_ii_svd$u, 2, est_trend)
      attributes(proj_ii_svd$u_detrended)$dimnames <- NULL
    } 
    if(lev_images){ proj_ii_svd$v <- proj_ii_svd$v[,chosen_PCs_ordered,drop=FALSE] }
    projection[[proj_ii_name]] = list(indices=chosen_PCs, svd=proj_ii_svd)
  }

  # ----------------------------------------------------------------------------
  # For each method... ---------------------------------------------------------
  # ----------------------------------------------------------------------------

  for(ii in 1:length(methods)){

    # --------------------------------------------------------------------------
    # Measure outlingness and ID outliers. -------------------------------------
    # --------------------------------------------------------------------------

    method_ii <- methods[ii]
    method_ii_split <- unlist(strsplit(method_ii, "__"))
    proj_ii_name <- method_ii_split[1]; out_ii_name <- method_ii_split[2]
    proj_ii <- projection[[proj_ii_name]]
    
    if (verbose) { 
      cat(paste0("Method ", method_ii)) 
      if (get_outliers) { cat(":") }
    }

    # Adjust PC number if using robust distance.
    if (out_ii_name %in% c("robdist", "robdist_bootstrap")) {
      # Let q <- N_/T_ (ncol(U)/nrow(U)). robustbase::covMcd()
      #  requires q <= approx. 1/2 for computation of the MCD covariance estimate.
      #  Higher q will use more components for estimation, thus retaining a
      #  higher resolution of information. Lower q will have higher breakdown
      #  points, thus being more resistant to outliers. (The maximal breakdown
      #  value is (N_ - T_ + 2)/2.) Here, we select q = 1/3 to yield a breakdown
      #  value of approx. 1/3. 
      #
      # For computational feasibility, we also limit the total number of PCs to 
      # 100.
      q <- 1/3
      max_keep = min(100, ceiling(T_*q))
      if (max_keep < length(proj_ii$indices)) {
        cat(paste(
          " Reducing number of PCs from", length(proj_ii$indices), "to", max_keep
        ))

        # Identify the indices to keep.
        if (proj_ii_name == "PCA_kurt") {
          max_idx <- sort(proj_ii$indices)[max_keep]
          proj_ii$indices <- proj_ii$indices[proj_ii$indices <= max_idx]
        } else {
          proj_ii$indices <- proj_ii$indices[1:max_keep]
        }

        # Take that SVD subset.
        proj_ii$svd$u <- proj_ii$svd$u[,1:max_keep,drop=FALSE]
        proj_ii$svd$d <- proj_ii$svd$d[1:max_keep]
        if ("PCs" %in% detrend && proj_ii_name!="PCATF") {
          proj_ii$svd$u_detrended <- proj_ii$svd$u_detrended[,1:max_keep,drop=FALSE]
        }
        if (solve_dirs) {
          proj_ii$svd$v <- proj_ii$svd$v[,1:max_keep,drop=FALSE]
        }
      }
    }

    # Compute the outlyingness measure.
    U_meas <- ifelse(
      "PCs" %in% detrend && proj_ii_name != "PCATF", 
      "u_detrended", 
      "u"
    )
    U_meas <- proj_ii$svd[[U_meas]]

    out_fun_ii <- switch(out_ii_name,
      leverage = out_measures.leverage,
      #robdist_bootstrap = out_measures.robdist_bootstrap,
      robdist = out_measures.robdist,
    )
    if (get_outliers) {
      cutoff_ii <- switch(out_ii_name,
        leverage=lev_cutoff,
        #robdist_bootstrap=rbd_cutoff,
        robdist=rbd_cutoff
      )
    } else {
      cutoff_ii <- NULL
    }
    out_kwargs_ii <- switch(out_ii_name,
      leverage = list(median_cutoff=cutoff_ii),
      #robdist_bootstrap = list(R_true=R_true, quantile_cutoff=cutoff_ii),
      robdist = list(quantile_cutoff=cutoff_ii)
    )
    out_kwargs_ii <- c(list(U = U_meas), out_kwargs_ii)
    out_ii <- do.call(out_fun_ii, out_kwargs_ii)
    outlier_measures[[method_ii]] <- out_ii$meas
    if (get_outliers) {
      outlier_cutoffs[[method_ii]] <- out_ii$cut
      outlier_flags[[method_ii]] <- out_ii$flag
    }
    if (out_ii_name == "robdist") {
      robdist_info[[method_ii]] <- out_ii$info
    }
      
    # --------------------------------------------------------------------------
    # Make leverage images.-----------------------------------------------------
    # --------------------------------------------------------------------------

    if (get_outliers) {
      if (sum(out_ii$flag) > 0) {
        if (verbose) {
          cat(" Outliers detected.")
          if (lev_images) {
            cat(" Computing leverage images.")
            outlier_lev_imgs[[method_ii]] <- get_leverage_images(
              proj_ii$svd, which(out_ii$flag), const_mask
            )
          }
        }
      } else {
        if (verbose) {
          cat(" No outliers detected.")
          if (lev_images) {
            cat(" Skipping leverage images.")
          }
        }
        outlier_lev_imgs[[method_ii]] <- list(mean=NULL, top=NULL, top_dir=NULL)
      }
    }

    if (verbose) {cat("\n")}
  }

  # ----------------------------------------------------------------------------
  # Format output. -------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # Conserve memory.
  gc()

  if (length(robdist_info) < 1) {robdist_info <- NULL}

  # Organize the output.
  if (verbose) { cat("Done!\n\n") }
  result <- list(
    params = NULL, 
    projections = projection, 
    outlier_measures = outlier_measures
  )
  if ("robdist" %in% out_meas) {
    result$robdist_info <- robdist_info
  }
  if (get_outliers) {
    result$outlier_cutoffs <- outlier_cutoffs
    result$outlier_flags <- outlier_flags
  }
  if (lev_images) { result$outlier_lev_imgs <- outlier_lev_imgs }
  result$params <- list(
    projection = names(projection),
    out_meas = out_meas,
    DVARS = DVARS,
    detrend = detrend,
    PCATF_kwargs = PCATF_kwargs,
    kurt_quantile = kurt_quantile,
    get_outliers = get_outliers,
    lev_cutoff = lev_cutoff,
    rbd_cutoff = rbd_cutoff,
    lev_images = lev_images,
    verbose = verbose
  )

  structure(result, class="clever")
}