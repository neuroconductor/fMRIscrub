#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' @param X A wide (observations x variables) numerical data matrix.
#' @param projection_methods Vector of projection methods to use. Choose at least
#'  one of the following: \code{"PCA_var"} for PCA + variance, \code{"PCA_kurt"}
#'  for PCA + kurtosis, and \code{"PCATF"} for PCA Trend Filtering + variance.
#'  Or, use \code{"all"} to use all methods. Default is \code{c("PCA_kurt")}.
#' 
#'  clever will use all combinations of the requested projection and 
#'  outlyingness methods that make sense. For example, if  
#'  \code{projection_methods=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{outlyingness_methods=c("leverage", "robdist")} then these five
#'  combinations will be used: PCATF with leverage, PCA_var with leverage, 
#'  PCA_var with robdist, PCA_kurt with leverage, and PCA_kurt with robdist. 
#'  Each method combination will yield its own outlyingness time series.
#' @param outlyingness_methods Vector of outlyingness measures to compute. Choose
#'  at least one of the following: \code{"leverage"} for leverage, or 
#'  \code{"robdist"} for robust distance. Or, use \code{"all"} to use all methods. 
#'  Default is \code{c("leverage")}.
#' 
#'  clever will use all combinations of the requested projection and 
#'  outlyingness methods that make sense. For example, if  
#'  \code{projection_methods=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{outlyingness_methods=c("leverage", "robdist")} then these five
#'  combinations will be used: PCATF with leverage, PCA_var with leverage, 
#'  PCA_var with robdist, PCA_kurt with leverage, and PCA_kurt with robdist. Each 
#'  method combination will yield its own outlyingness time series.
#' @param DVARS Should DVARS (Afyouni and Nichols, 2017) be computed too? Default 
#'  is \code{TRUE}.
#' @param PCs_detrend Detrend all PCs before computing leverage or robust distance? 
#'  Detrending is recommended for time-series data, especially if there are many
#'  time points or changing circumstances, such as in task-based fMRI. Detrending
#'  should not be used with non-time-series data because the observations are 
#'  not temporally related. Default is \code{TRUE}.
#' @param PCATF_kwargs Named list of arguments for PCATF: maximum number of PCs
#'  to compute \code{K} (Default 1000) the trend filtering parameter lambda
#'  (Default \code{0.5}), the number of iterations \code{niter_max} (Default 
#'  \code{1000}), convergence tolerance \code{tol} (Default \code{1e-8}), and option 
#'  to print updates \code{verbose} (Default \code{FALSE}).
#' @param kurt_quantile What cutoff quantile for kurtosis should be used? This
#'  argument only applies if \code{"PCA_kurt"} is one of the projection methods 
#'  in \code{"projection_methods"}. Default is \code{0.95}.
#' @param kurt_detrend Should the PCs be detrended before measuring kurtosis? 
#'  This argument only applies if \code{"PCA_kurt"} is one of the projection methods 
#'  in \code{"projection_methods"}. Default is \code{TRUE}. Detrending is highly
#'  recommended for time-series data, because trends can induce high kurtosis even
#'  in the absence of outliers. Detrending should not be done with non-time-series
#'  data because the observations are not temporally related.
#' @param id_outliers Should the outliers be identified? Default is \code{TRUE}.
#' @param lev_cutoff The outlier cutoff value for leverage, as a multiple of the median
#'  leverage. Only used if 
#'  \code{"leverage" \%in\% projection_methods} and \code{id_outliers}. Default is 
#'  \code{4}, or \eqn{4 * median}.
#' @param rbd_cutoff The outlier cutoff quantile for MCD distance. Only used if 
#'  \code{"robdist" \%in\% projection_methods} and \code{id_outliers}. Default 
#"  is \code{0.9999}, for the \eqn{0.9999} quantile.
#'  The quantile is computed from the estimated F distribution.
#' @param lev_images Should leverage images be computed? If \code{FALSE} memory is
#'  conserved. Default is \code{FALSE}.
#' @param verbose Should occasional updates be printed? Default is \code{FALSE}.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#'  \item{params}{A list of all the arguments used.}
#'  \item{projections}{
#'    \describe{
#'      \item{PC_var}{
#'        \describe{
#'          \item{indices}{The indices retained from the original SVD decomposition 
#'            to make the variance-based PC projection.} 
#'          \item{PCs}{The subsetted SVD decomposition.}  
#'        }
#'      }
#'      \item{PC_kurt}{
#'        \describe{
#'          \item{indices}{The indices retained from the original SVD decomposition 
#'            to make the kurtosis-based PC projection. They are ordered from highest 
#'            kurtosis to lowest kurtosis.}  
#'          \item{PCs}{The subsetted SVD decomposition. PCs are ordered in the standard
#'            way, from highest variance to lowest variance, instead of by kurtosis.}  
#'        }
#'      }
#'      \item{PCATF}{
#'        \describe{
#'          \item{indices}{The indices of the trend-filtered PCs used to make the
#'            projection.}  
#'          \item{PCs}{The PCATF result.}  
#'        }
#'      }
#'    }
#'  }
#'  \item{outlier_measures}{
#'    \describe{
#'      \item{PC_var__lev}{The leverage values for the PC_var projection.}
#'      \item{PC_kurt__lev}{The leverage values for the PC_kurt projection.}
#'      \item{PCATF__lev}{The leverage values for the PCATF projection.}
#'      \item{PC_var__rbd}{The robust MCD distance values for the PC_var projection.}
#'      \item{PC_kurt__rbd}{The robust MCD distance values for the PC_kurt projection.}
#'      \item{DVARS_DPD}{The Delta percent DVARS values.}
#'      \item{DVARS_ZD}{The DVARS z-scores.}
#'    }
#'  }
#'  \item{outlier_cutoffs}{
#'    \describe{
#'      \item{lev}{The leverage cutoff for outlier detection: \code{lev_cutoff} times
#'        the median leverage.}
#'      \item{MCD}{The robust distance cutoff for outlier detection: the 
#'        \code{rbd_cutoff} quantile of the estimated F distribution.}
#'      \item{DVARS_DPD}{The Delta percent DVARS cutoff: +/- 5 percent}
#'      \item{DVARS_ZD}{The DVARS z-score cutoff: the one-sided 5 percent 
#'        significance level with Bonferroni FWER correction.}
#'    }
#'  }
#'  \item{outlier_flags}{
#'    \describe{
#'      \item{PC_var__leverage}{Whether each observation surpasses the outlier cutoff.}
#'      \item{PC_kurt__leverage}{Whether each observation surpasses the outlier cutoff.}
#'      \item{PCATF__leverage}{Whether each observation surpasses the outlier cutoff.}
#'      \item{PC_var__robdist}{Whether each observation surpasses the outlier cutoff.}
#'      \item{PC_kurt__robdist}{Whether each observation surpasses the outlier cutoff.}
#'      \item{DVARS_DPD}{Whether each observation surpasses the outlier cutoff.}
#'      \item{DVARS_ZD}{Whether each observation surpasses the outlier cutoff.}
#'    }
#'  }
#'  \item{MCD_scale}{The scale value for out-of-MCD observations, and NA for
#'    in-MCD observations. NULL if \code{method} is not robust distance.}
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
  projection_methods = "PCA_kurt",
  outlyingness_methods = "leverage",
  DVARS = TRUE,
  PCs_detrend = TRUE,
  PCATF_kwargs = NULL,
  kurt_quantile = .95,
  kurt_detrend = TRUE,
  id_outliers = TRUE,
  lev_cutoff = 4,
  rbd_cutoff = 0.9999,
  lev_images = FALSE,
  verbose = FALSE) {

  ###################
  # PRELIMINARY STEPS
  ###################

  # Define the projection and outlyingness methods, and their valid combos.
  all_projection_methods <- c("PCA_var", "PCA_kurt", "PCATF")
  all_outlyingness_methods <- c("leverage", "robdist")
  all_valid_methods <- c(
    c("PCA_var__leverage", "PCA_kurt__leverage", 
      "PCA_var__robdist", "PCA_kurt__robdist", 
      "PCATF__leverage")
  )

  TOL <- 1e-8 # cutoff for detection of zero variance/MAD voxels

  # Check arguments.
  if(!is.matrix(X)){ X <- as.matrix(X) }
  if(identical(projection_methods, "all")){
    projection_methods <- all_projection_methods
  } else {
    projection_methods <- match.arg(
      projection_methods, all_projection_methods, 
      several.ok=TRUE)
  }
  if(identical(outlyingness_methods, "all")){
    outlyingness_methods <- all_outlyingness_methods
  } else {
    outlyingness_methods <- match.arg(
      outlyingness_methods, all_outlyingness_methods, 
      several.ok=TRUE)
  }
  stopifnot(is.logical(DVARS))
  stopifnot(is.logical(PCs_detrend))
  if(!identical(PCATF_kwargs, NULL)){
    names(PCATF_kwargs) <- match.arg(
      names(PCATF_kwargs), c("K", "lambda", "niter_max", "tol", "verbose"),
      several.ok=TRUE)
    if(length(PCATF_kwargs) != length(unique(unlist(PCATF_kwargs)))){
      stop("Duplicate PCATF_kwargs were given.\n")
    }
  }
  if("PCA_kurt" %in% projection_methods){
    stopifnot(is.numeric(kurt_quantile))
    stopifnot((kurt_quantile < 1) & (kurt_quantile > 0))
    stopifnot(is.logical(kurt_detrend))
  }
  stopifnot(is.logical(id_outliers))
  if((lev_images) & (!id_outliers)){
    stop("Invalid argument: computing leverage images requires id_outliers==TRUE.\n")
  }
  stopifnot(is.numeric(lev_cutoff))
  stopifnot(lev_cutoff > 0)
  stopifnot((rbd_cutoff > 0) & (rbd_cutoff < 1))
  stopifnot(is.logical(lev_images))
  stopifnot(is.logical(verbose))

  Npre_ <- ncol(X)
  T_ <- nrow(X)
  if(Npre_ < T_){
    warning("Data matrix has more rows than columns. Check that observations\
             are in rows and variables are in columns.\n")
  }

  # Collect all the methods to compute.
  methods <- all_valid_methods[all_valid_methods %in% 
    outer(projection_methods, outlyingness_methods, paste, sep='__')
  ]
  if(length(methods) < 1){
    stop("No valid method combinations. Check that the projection and\
          outlyingness methods are compatible.\n")
  }
  outlier_measures <- outlier_lev_imgs <- setNames(vector("list", length(methods)), methods)
  if(id_outliers){
    outlier_cutoffs <- outlier_flags <- setNames(vector("list", length(methods)), methods)
  }
  MCD_scale <- vector("list")

  # Center and scale the data robustly.
  # Do it here instead of calling scale_med to save memory.
  if(verbose){ print("Centering and scaling the data matrix.") }
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
      " ). These will be removed for estimation of the covariance.\n"))
    }
  }
  mad <- mad[!const_mask]
  X <- X[!const_mask,]
  X <- X/c(mad)
  # Revert transpose.
  X <- t(X)
  N_ <- ncol(X)
   
  # Compute DVARS.
  if(DVARS){
    if(verbose){ print("Computing DVARS.") }
    X_DVARS <- compute_DVARS(X, normalize=FALSE, norm_I=100, verbose=verbose)
    
    outlier_measures$DVARS_DPD <- X_DVARS$DPD
    outlier_measures$DVARS_ZD <- X_DVARS$ZD

    if(id_outliers){
      outlier_cutoffs$DVARS_DPD <- 5
      outlier_cutoffs$DVARS_ZD <- qnorm(1-.05/T_)
      
      outlier_flags$DVARS_DPD <- X_DVARS$DPD > outlier_cutoffs$DVARS_DPD
      outlier_flags$DVARS_ZD <- X_DVARS$ZD > outlier_cutoffs$DVARS_ZD
      outlier_flags$DVARS_both <- outlier_flags$DVARS_DPD & outlier_flags$DVARS_ZD
    }
  }

  ##################
  # DATA PROJECTIONS
  ##################

  # Compute the PC scores (and directions, if leverage images or PCATF are desired).
  solve_directions <- (lev_images) | ("PCATF" %in% projection_methods)
  if(verbose){
    print(paste0(
      "Computing the",
      ifelse(("PCA_var" %in% projection_methods) | ("PCA_kurt" %in% projection_methods),
        ifelse("PCATF" %in% projection_methods, " normal and trend-filtered", ""),
        ifelse("PCATF" %in% projection_methods, " trend-filtered", "INTERNAL ERROR")),
      " PC scores",
      ifelse(solve_directions, " and directions", ""), "."))
  }
  if(solve_directions){
    X.svd <- svd(X)
    if(!("PCATF" %in% projection_methods)){ rm(X) }
  } else {
    # Avoid computing U matrix to conserve memory.
    XXt <- tcrossprod(X)
    if(!("PCATF" %in% projection_methods)){ rm(X) }
    X.svd <- svd(XXt)
    rm(XXt)
    X.svd$d <- sqrt(X.svd$d)
    X.svd$v <- NULL
  }

  # Compute PCATF, if requested.
  if("PCATF" %in% projection_methods){
    X.svdtf <- do.call(PCATF, 
      c(list(X=X, X.svd=X.svd, solve_directions=solve_directions), PCATF_kwargs))
    # The PC directions were needed to compute PCATF. If leverage images 
    #   are not wanted, we can now delete the directions to save space.
    if(!lev_images){ X.svd$v <- NULL }

    # Remove trend-filtered PCs with constant scores.
    tf_zero_var <- apply(X.svdtf$u, 2, var) < TOL
    if(any(tf_zero_var)){
      if(all(tf_zero_var)){
        stop("Error: All trend-filtered PC scores are zero-variance.\n")
      }
      warning(paste("Warning:", sum(tf_zero_var), 
        "trend-filtered PC scores are zero-variance. Removing these PCs.\n"))
      X.svdtf$u <- X.svdtf$u[,!tf_zero_var]
      X.svdtf$d <- X.svdtf$d[!tf_zero_var]
      if(lev_images){ X.svdtf$v <- X.svdtf$v[,!tf_zero_var] }
    }
  }
  gc()

  # Choose which PCs to retain for each projection.
  projections <- setNames(vector("list", length(projection_methods)), projection_methods)
  for(i in 1:length(projection_methods)){
    projection_method <- projection_methods[i]

    if(verbose){
      print(switch(projection_method,
        PCA_var = "Identifying the PCs with high varaince.",
        PCA_kurt = "Identifying the PCs with high kurtosis.",
        PCATF = "Identifying the trend-filtered PCs with high varaince."
      ))
    }

    if(projection_method=="PCATF"){
      # We have already computed the PCs we want.
      if("max_keep" %in% names(choose_PCs_kwargs)){
        chosen_PCs <- 1:min(ncol(X.svdtf$u), choose_PCs_kwargs$max_keep)
      } else {
        chosen_PCs <- 1:ncol(X.svdtf$u)
      }
    } else {
      choose_PCs_kwargs <- list(svd=X.svd)
      if(projection_method == "PCA_var"){
        chosen_PCs = choose_PCs.variance(svd=X.svd)
      }else if(projection_method == "PCA_kurt"){
        chosen_PCs = choose_PCs.kurtosis(
          svd=X.svd, kurt_quantile=kurt_quantile, detrend=kurt_detrend)
      }
    }
    chosen_PCs_ordered <- chosen_PCs[order(chosen_PCs)] # kurtosis order =/= index order
    
    projection = list(indices = chosen_PCs)
    if(projection_method=="PCATF"){
      projection$svd <- X.svdtf
    } else {
      projection$svd <- X.svd
    }
    projection$svd$u <- projection$svd$u[,chosen_PCs_ordered]
    if(PCs_detrend & projection_method!="PCATF"){
      projection$svd$u_detrended <- projection$svd$u - apply(projection$svd$u, 2, est_trend)
      attributes(projection$svd$u_detrended)$dimnames <- NULL
    } 
    projection$svd$d <- projection$svd$d[chosen_PCs_ordered]
    if(lev_images){ projection$svd$v <- projection$svd$v[,chosen_PCs_ordered] }
    projections[[projection_method]] = projection
  }

  ##########################################
  # OUTLIER MEASUREMENTS AND LEVERAGE IMAGES
  ##########################################

  # Perform the rest of the steps for each method combination.
  for(i in 1:length(methods)){
    method_combo <- methods[i]
    method_split <- unlist(strsplit(method_combo, "__"))
    projection_method <- method_split[1]
    outlyingness_method <- method_split[2]
    projection <- projections[[projection_method]]
    
    if(verbose){
      print(paste0("Method ", method_combo, ":"))
    }

    # Adjust PC number if using robust distance.
    if(outlyingness_method == "robdist"){
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
      if(max_keep < length(projection$indices)){
        print(paste("...reducing number of PCs from", length(projection$indices),
          "to", max_keep))
        if(projection_method == "PCA_kurt"){
          max_idx <- sort(projection$indices)[max_keep]
          projection$indices <- projection$indices[projection$indices <= max_idx]
        } else {
          projection$indices <- projection$indices[1:max_keep]
        }
        projection$svd$u <- projection$svd$u[,1:max_keep]
        if(PCs_detrend & projection_method!="PCATF"){
          projection$svd$u_detrended <- projection$svd$u_detrended[,1:max_keep]
        }
        projection$svd$d <- projection$svd$d[1:max_keep]
        if(solve_directions){
          projection$svd$v <- projection$svd$v[,1:max_keep]
        }
      }
    }

    # Compute outlyingness measure.
    outlyingness_method.fun <- switch(outlyingness_method,
      leverage = PC.leverage,
      robdist = PC.robdist
    )
    if(PCs_detrend & projection_method!="PCATF"){
      U_meas <- projection$svd$u_detrended
    } else {
      U_meas <- projection$svd$u
    }
    if(outlyingness_method == "robdist"){
      measure <- outlyingness_method.fun(U_meas)
      this_MCD_scale <- ifelse(measure$inMCD, NA, measure$outMCD_scale)
      Fparam <- measure$Fparam
      measure <- measure$mah
      MCD_scale[[method_combo]] <- this_MCD_scale
    } else {
      measure <- outlyingness_method.fun(U_meas)
    }
    outlier_measures[[method_combo]] <- measure
    
    # Identify outliers.
    if(id_outliers){
      if(verbose){ print("...identifying outliers.") }
      if(outlyingness_method == "leverage"){
        out_cutoff <- lev_cutoff * median(measure)
        out_flag <- measure > out_cutoff
      } else if(outlyingness_method == "robdist"){
        out_cutoff <- qf(p=rbd_cutoff, df1=Fparam$df[1], df2=Fparam$df[2])
        out_flag <- measure*this_MCD_scale > out_cutoff
        out_flag[is.na(out_flag)] <- FALSE
      } else {
        stop("INTERNAL ERROR: outlyingness_method not recognized.")
      }
      outlier_cutoffs[[method_combo]] <- out_cutoff
      outlier_flags[[method_combo]] <- out_flag
    }

    # Make leverage images.
    if(lev_images){
      if(sum(out_flag) > 0){
        if(verbose){
          print("...outliers detected. Computing leverage images.")
        }
        outlier_lev_imgs[[method_combo]] <- get_leverage_images(
          projection$svd, which(out_flag), const_mask)
      } else {
        if(verbose){
          print("...no outliers detected. Skipping leverage images.")
        }
        outlier_lev_imgs[[method_combo]] <- list(mean=NULL, top=NULL, top_dir=NULL)
      }
    }
  }

  ###############
  # FORMAT OUTPUT
  ###############

  # Conserve memory.
  gc()

  if(length(MCD_scale) < 1){MCD_scale <- NULL}

  # Organize the output.
  if(verbose){ print("Done! Organizing results.") }
  result <- list(
    params = NULL, 
    projections = projections, 
    outlier_measures = outlier_measures)
  if("robdist" %in% outlyingness_methods){
    result$MCD_scale <- MCD_scale
  }
  if(id_outliers){
    result$outlier_cutoffs <- outlier_cutoffs
    result$outlier_flags <- outlier_flags
  }
  if(lev_images){ result$outlier_lev_imgs <- outlier_lev_imgs }
  result$params <- list(
    projection_methods = projection_methods,
    outlyingness_methods = outlyingness_methods,
    DVARS = DVARS,
    PCs_detrend = PCs_detrend,
    PCATF_kwargs = PCATF_kwargs,
    kurt_quantile = kurt_quantile,
    kurt_detrend = kurt_detrend,
    id_outliers = id_outliers,
    lev_cutoff = lev_cutoff,
    rbd_cutoff = rbd_cutoff,
    lev_images = lev_images,
    verbose = verbose)

  class(result) <- c("clever", class(result))

  return(result)
}