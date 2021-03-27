#' Order of operations in \code{clever}
#' 
#' @section Order of operations:
#'  The input data undergo a nuisance regression to remove low frequency trends
#'  (see the \code{DCT} argument to control or disable this behavior) and any
#'  other nuisance signals (see the \code{nuisance_too} argument to specify them).
#'  (If \code{DCT==0} and \code{is.null(nuisance_too)}, no nuisance regression 
#'  is performed.) Then the data are centered and scaled (see the \code{center}
#'  and \code{scale} arguments to control or disable this behavior). Leverage
#'  is computed afterward.
#'  
#' @name clever_order_of_operations
#' @keywords internal
NULL

#' Removed arguments
#' 
#' @param measures Character vector indicating the measures to compute. Choose 
#'  at least one of the following: 
#' 
#'  \describe{
#'    \item{\code{"leverage"}}{Leverage scrubbing, which is based on
#'      projecting the data onto directions thought to express outlier 
#'      information.}
#'    \item{\code{"robdist"}}{Robust Mahalanobis-based distance}
#'    \item{\code{"DVARS"}}{Traditional DVARS}
#'    \item{\code{"DVARS2"}}{Delta-percent-DVARS and z-score-DVARS (Afyouni and 
#'      Nichols, 2018)}
#'    \item{\code{"FD"}}{Framewise Displacement. Requires \code{motion}.}
#'    \item{\code{"motion"}}{Translation and rotation realignment parameters. 
#'      Requires \code{motion}.}
#'    \item{\code{"GSR"}}{Global Signal of the data.}
#'  }
#' 
#'  Use \code{"all"} to select all available measures. (FD and motion will only 
#'  be computed if the motion realignment parameters are provided.) Default: 
#'  \code{"leverage", "DVARS2"}.
#'
#'  Note that motion and GSR are not direct measures of outlyingness,
#'  so they do not have corresponding \code{outlier_cutoffs}.
#' @param motion Only used if the \code{"FD"} measure is requested. An 
#'  \eqn{N \times 6} matrix in which the first three columns represent the
#'  translational realignment parameters (mm), and the second three columns represent
#'  the rotational realignment parameters in (radians). To convert radians to mm,
#'  the displacement on a sphere of radius 50 mm will be computed.
#'
#'  Alternatively, this can be the file path to an \eqn{N \times 6} matrix which can be
#'  read with \code{\link{read.table}} (fields separated by white-space; no
#'  header).
#' @param noise_nPC Only applies to CompCor.
#'  The number of principal components to compute for each noise
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
#' @param noise_erosion Only applies to CompCor.
#'  The number of voxel layers to erode the noise ROIs by. 
#' 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its 
#'  elements will be recycled. Default: \code{NULL}, which will use a value of
#'  0 (do not erode the noise ROIs).
#' @param DCT Detrend the columns of the data using the discrete cosine
#'  transform (DCT)? Use an integer to indicate the number of cosine bases to 
#'  use for detrending. Use \code{0} (default) to forgo detrending. 
#' 
#'  The data must be centered, either before input or with \code{center}.
#' 
#'  Detrending is highly recommended for time-series data, especially if there 
#'  are many time points or evolving circumstances affecting the data. Additionally,
#'  if kurtosis is being used to select the projection directions, trends can 
#'  induce positive or negative kurtosis, contaminating the connection between 
#'  high kurtosis and outlier presence. 
#'  
#'  Detrending should not be used with non-time-series data because the 
#'  observations are not temporally related.
#' @param PESEL Leverage is based on the largest \eqn{k} PCA, ICA, or PCATF 
#'  components. \code{k} can be determined by PESEL with \code{\link[pesel]{pesel}}
#'  (default, \code{PESEL==TRUE}). Otherwise, \eqn{k} will be the number of 
#'  above-average-variance PCs (\code{PESEL==FALSE}). In either case, \code{k}
#'  is calculated based on PCA, even to select the number of ICA and PCATF
#'  components.  
#' 
#'  Note that not all \code{k} components are ultimately used: only components with
#'  high kurtosis (> 99th percentile of the kurtosis of data of equal length 
#'  from a Normal distribution) contribute to the leverage calculation.
#' @param nuisance_too A matrix of nuisance signals to regress from the data
#'  before, i.e. a "design matrix." Should have \eqn{T} rows. Nuisance
#'  regression will be performed simultaneously with DCT detrending if 
#'  applicable. \code{NULL} (default) to not add additional nuisance regressors.
#' @name removed_arguments
#' @keywords internal
NULL

#' clever
#' 
#' @param X Wide numeric data matrix (\eqn{T observations \times V variables}, 
#'  \eqn{T << V}). If \code{X} represents an fMRI run, \eqn{T} should be the 
#'  number of timepoints and \eqn{V} should be the number of brainordinate 
#'  vertices/voxels. The outlyingness of the rows in \code{X} will be measured.
#' @param nuisance Nuisance signals to regress from each column of \code{X}.
#'  Should be a \eqn{T \times N} numeric matrix where \eqn{N} represents the
#'  number of nuisance signals. Default: a matrix with a constant column
#'  (this is the intercept term in the design matrix) and four DCT bases. 
#'  This default nuisance regression will have the effect of detrending the
#'  data by removing low-frequency components. To not perform any nuisance
#'  regression set this argument to \code{NULL}, \code{0}, or \code{FALSE}.
#' 
#'  Detrending is highly recommended for time-series data, especially if there 
#'  are many time points or evolving circumstances affecting the data. Additionally,
#'  if kurtosis is being used to select the projection directions, trends can 
#'  induce positive or negative kurtosis, contaminating the connection between 
#'  high kurtosis and outlier presence. Detrending should not be used with 
#'  non-time-series data because the observations are not temporally related.
#' 
#'  To perform nuisance regression with an intercept, DCT bases, and additional
#'  signals, specify something like \code{cbind(1, dct_bases(nrow(X), 4), more_nuisance)}.
#' @param center,scale Center the columns of the data by their medians, and scale the
#'  columns of the data by their median absolute distances (MADs)? Default: \code{TRUE}. 
#'  Centering is necessary for computing the projections, so if \code{center} is
#'  \code{FALSE}, the data must already be centered.
#' 
#'  Note that centering and scaling occur after nuisance regression, so even if
#'  \code{center} is \code{FALSE}, the data will be centered on the means if
#'  the nuisance regression included an intercept term, as it does by default.
#' @param var_detrend Stabilize the variance of the PCA, PCATF, and ICA components
#'  prior to computing leverage? \code{TRUE} (default), \code{FALSE}, or the number
#'  of DCT bases to use for variance detrending (\code{TRUE} will use 2). 
#'  Note that variance stabilization is performed on the projection scores and 
#'  not the data itself.
#' 
#'  For each component, a linear model of four DCT bases and an intercept is fit
#'  on the log squared values after demeaning. The score timeseries of this 
#'  component is divided by the smooth, estimated variance timeseries from this 
#'  linear model. Finally, the values are re-centered and scaled to match the original
#'  mean and variance, yielding the variance stabilized scores. 
#' @param kurt_quantile Only applies to PCA and ICA leverage. What cutoff quantile
#'  for kurtosis should be used to select the components? Default: \code{0.99}.
#' @param PCATF_kwargs Arguments to \code{\link{PCATF}} in list form. Valid
#'  entries are:
#'  
#'  \describe{
#'    \item{lambda}{Trend-filtering parameter. Default: \code{5}.}
#'    \item{niter_max}{Maximum number of iterations. Default: \code{1000}.}
#'    \item{TOL}{Convergence tolerance parameter. Default: \code{1e-8}.}
#'    \item{verbose}{Print updates? Default: \code{FALSE}.}
#'  }
#' @param get_dirs Do the projection directions need to be returned? This is the 
#'  \eqn{V} matrix in PCA and \eqn{S} matrix in ICA. The default is \code{FALSE}
#'  to save memory. However, \code{get_dirs==TRUE} is required for \code{\link{lev_images}}.
#' @param full_PCA Only applies to the PCA projection. Return the full SVD? 
#'  Default: \code{FALSE} (return only the components used to compute the measures).
#' @param get_outliers Should outliers be flagged based on \code{cutoff}? Default: \code{TRUE}.
#' @param cutoff Median leverage cutoff value. Default: \code{4}.
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @name clever_Params
#' @keywords internal
NULL

#' fMRI data for clever and CompCor
#' 
#' @param X Wide numeric data matrix (\eqn{T observations \times V variables}, \eqn{T << V}).
#'  For example, if \code{X} represents an fMRI run, \eqn{T} should be the number
#'  of timepoints and \eqn{V} should be the number of brainordinate vertices/voxels.
#' 
#'  Or, a 4D array or NIFTI or file path to a NIFTI (\eqn{I \times J \times K \times T} 
#'  observations), in which case \code{ROI_data} must be provided. 
#'  (The vectorized data will be \eqn{T timepoints \times V_{in-mask} voxels})
#' 
#'  Or, a \code{ciftiTools} code{"xifti"} object or a file path to a CIFTI
#'  (The vectorized data will be \eqn{T timepoints \times V_{left+right+sub} greyordinates}).
#' @param ROI_data Indicates the data ROI. Allowed arguments depend on \code{X}:
#' 
#'  If \code{X} is a matrix, this must be a length \eqn{V} logical vector, where
#'  the data ROI is indicated by \code{TRUE} values. If \code{"infer"} (default), all 
#'  columns of \code{X} will be included in the data ROI (\code{rep(TRUE, V)}).
#' 
#'  If \code{X} is an array or NIFTI, this must be either a vector of values
#'  to expect for out-of-mask voxels in \code{X}, or a (file path to a) 3D NIFTI.
#'  In the latter case, each of the volume dimensions should match the first
#'  three dimensions of \code{X}. Voxels in the data ROI should be indicated by
#'  \code{TRUE} and all other voxels by \code{FALSE}. If \code{"infer"} (default),
#'  will be set to \code{c(0, NA, NaN)} (include all voxels which are not constant
#'  \code{0}, \code{NA}, or \code{NaN}).
#' 
#'  If \code{X} is a \code{"xifti"} this must be the \code{brainstructures}
#'  argument to \code{\link[ciftiTools]{read_xifti}}. If \code{"infer"} (default),
#'  \code{brainstructures} will be set to \code{"all"} (use both left and right
#'  cortex vertices, and subcortical voxels).
#'
#'  If \code{NULL}, the data ROI will be empty. This is useful for obtaining just
#'  the noise ROI, if the data and noise are located in separate files.
#' @param ROI_noise Indicates the noise ROIs for aCompCor. Should be a list where
#'  each entry corresponds to a distinct noise ROI. The names of the list should
#'  be the ROI names, e.g. \code{"white_matter"} and \code{"csf"}. The expected
#'  formats of the list entries depends on \code{X}:
#' 
#'  For all types of \code{X}, \code{ROI_noise} entries can be a matrix of noise
#'  ROI data. The matrix should have \eqn{T} rows, with each column being a
#'  data location's timeseries.
#'  
#'  If \code{X} is a matrix, entries can also indicate a noise ROI within \code{X}.
#'  These entries must be a length \eqn{V} logical vector with \code{TRUE} values 
#'  indicating locations in \code{X} within that noise ROI. Since the ROIs must 
#'  not overlap, the masks must be mutually exclusive with each other, and with 
#'  \code{ROI_data}. 
#' 
#'  If \code{X} is an array or NIFTI, entries can also indicate a noise ROI within \code{X}.
#'  These entries must be a logical array or (file path to) a 3D NIFTI with the
#'  same spatial dimensions as \code{X}, and with \code{TRUE} values indicating
#'  voxels inside the noise ROI. Since the ROIs must not overlap, the masks must
#'  be mutually exclusive with each other, and with \code{ROI_data}. 
#' 
#'  (If \code{X} is a \code{"xifti"}, entries must be data matrices, since no 
#'  greyordinate locations in \code{X} are appropriate noise ROIs).
#' @name data_clever_CompCor_Params
NULL

#' noise parameters for CompCor
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
#' @name noise_Params
NULL