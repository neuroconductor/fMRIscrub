#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' @param X A wide (observations x variables) numerical data matrix.
#' @param projection_methods Vector of projection methods to use. Choose at least
#'  one of the following: \code{"PCA_var"} for PCA + variance, \code{"PCA_kurt"}
#'  for PCA + kurtosis, and \code{"PCATF"} for PCA Trend Filtering + variance.
#'  Default is \code{c("PCA_var")}.
#' 
#'  clever will use all combinations of the requested projection and 
#'  outlyingness methods that make sense. For example, if  
#'  \code{projection_methods=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{outlyingness_methods=c("leverage", "robdist")} then these four
#'  combinations will be used: PCATF with leverage, PCA_var with leverage, 
#'  PCA_var with robdist, and PCA_kurt with leverage. Each method combination 
#'  will yield its own outlyingness time series.
#' @param outlyingness_methods Vector of outlyingness measures to compute. Choose
#'  at least one of the following: \code{"leverage"} for leverage, 
#'  \code{"robdist"} for robust distance, and \code{"robdist_subset"} for robust
#'  distance subset.  Default is \code{c("leverage")}.
#' 
#'  clever will use all combinations of the requested projection and 
#'  outlyingness methods that make sense. For example, if  
#'  \code{projection_methods=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{outlyingness_methods=c("leverage", "robdist")} then these four
#'  combinations will be used: PCATF with leverage, PCA_var with leverage, 
#'  PCA_var with robdist, and PCA_kurt with leverage. Each method combination 
#'  will yield its own outlyingness time series.
#' @param DVARS Should DVARS (Afyouni and Nichols, 2017) be computed too? Default 
#'  is \code{TRUE}.
#' @param PCATF_kwargs Named list of arguments for PCATF: the trend filtering 
#'  parameter lambda (Default \code{0.5}), the number of iterations niter_max 
#'  (Default \code{1000}), convergence tolerance tol (Default \code{1e-8}), and 
#'  option to print updates verbose (Default \code{FALSE}).
#' @param kurt_quantile What cutoff quantile for kurtosis should be used? This
#'  argument only applies if \code{"PCA_kurt"} is one of the projection methods 
#'  in \code{"projection_methods"}. Default is \code{0.9}.
#' @param kurt_detrend Should the PCs be detrended before measuring kurtosis? 
#'  This argument only applies if \code{"PCA_kurt"} is one of the projection methods 
#'  in \code{"projection_methods"}. Default is \code{TRUE}. Detrending is highly
#'  recommended for time-series data, because trends can induce high kurtosis even
#'  in the absence of outliers. It is highly advised against for non-time-series
#'  data because the observations are not temporally related.
#' @param id_out Should the outliers be identified? Default is \code{TRUE}.
#' @param lev_img_level An integer between \code{0} and \code{3}. We can use the
#'  selected PCs to visualize artifact signals ("leverage images") at each
#'  outlying time point. Options \code{1} to \code{3} will return the leverage
#'  images for time points meeting each respective threshold, with \code{1}
#'  being the lowest and \code{3} being the strictest threshold. They require
#'  \code{id_out} to be \code{TRUE}. Option \code{0} will not compute the
#'  leverage images, conserving memory. Option \code{1} is the default.
#'  \code{FALSE} will yield Option \code{0} and \code{TRUE} will yield Option
#'  \code{1}.
#' @param verbose Should occasional updates be printed? Default is \code{FALSE}.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#'   \item{params}{A list of all the arguments used.}
#'  \item{PCs}{
#'    \describe{
#'      \item{indices}{The original indices of the selected PCs.}
#'      \item{svd}{The selected subset of the SVD. The v matrix (PC directions)
#'        is witheld to conserve memory.}
#'    }
#'  }
#'  \item{leverage}{The leverage of each observation. NULL if \code{method} is
#'    not PCA leverage.}
#'  \item{robdist}{The robust distance of each observation. NULL if
#'    \code{method} is not robust distance (subset).}
#'  \item{inMCD}{Whether each observation is within the MCD subset. NULL if
#'    \code{method} is not robust distance (subset).}
#'  \item{outliers}{An n X 3 data.frame indicating if each observation is an
#'    outlier at each of the three levels.}
#'  \item{cutoffs}{Outlier cutoff values.}
#'  \item{lev_imgs}{
#'    \describe{
#'      \item{mean}{The average of the PC directions, weighted by the unscaled
#'        PC scores at each outlying time point (U[i,] * V^T). Row names are
#'        the corresponding time points.'}
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
  projection_methods = "PCA_var",
  outlyingness_methods = "leverage",
  DVARS = TRUE,
  PCATF_kwargs = NULL,
  kurt_quantile = .9,
  kurt_detrend = TRUE,
  id_out = TRUE,
  lev_img_level = 1,
  verbose = FALSE) {

  ###################
  # PRELIMINARY STEPS
  ###################

  # Define the projection and outlyingness methods, and their valid combos.
  all_projection_methods = c("PCA_var", "PCA_kurt", "PCATF")
  all_outlyingness_methods = c("leverage", "robdist", "robdist_subset")
  all_valid_method_combos = c(
    as.character(outer('PCA_var', all_outlyingness_methods, FUN=paste, sep='__')),
    c("PCA_kurt", "leverage"), c("PCATF", "leverage")
  )

  TOL <- 1e-8 # cutoff for detection of zero variance/MAD voxels

  # Check arguments.
  projection_methods <- match.arg(
    projection_methods,
    all_projection_methods, 
    several.ok=TRUE)
  outlyingness_methods <- match.arg(
    outlyingness_methods,
    all_outlyingness_methods, 
    several.ok=TRUE)
  if(!is.matrix(X)){ X <- as.matrix(X) }
  if("PCA_kurt" %in% projection_methods){
    if(!is.numeric(kurt_quantile)){
      stop("Invalid argument: kurt_quantile must be numeric.\n")
    }
    if((kurt_quantile > 1) | (kurt_quantile < 0)){
      stop("Invalid argument: kurt_quantile must be between 0 and 1.\n")
    }
    if(!is.logical(kurt_detrend)){
      stop("Invalid argument: kurt_detrend must be TRUE or FALSE.\n")
    }
  }
  if(!is.logical(id_out)){
    stop("Invalid argument: id_out must be TRUE or FALSE.\n")
  }
  if(!is.numeric(lev_img_level)){ lev_img_level <- as.numeric(lev_img_level) }
  if(!(lev_img_level %in% c(0,1,2,3))){
    stop("Invalid argument: lev_img_level must be 0, 1, 2, or 3.\n")
  }
  if((lev_img_level > 0) & (!id_out)){
    stop("Invalid argument: computing leverage images requires id_out==TRUE.\n")
  }
  if(!is.logical(verbose)){
    stop("Invalid argument: verbose must be TRUE or FALSE.\n")
  }

  method_combos = all_valid_method_combos[all_valid_method_combos %in% 
    outer(all_projection_methods, all_outlyingness_methods, paste, '__')
  ]

  projection.fun <- switch(projection_method, 
    kurtosis=choose_PCs.kurtosis,
    variance=choose_PCs.variance)
  outlyingness.fun <- switch(outlyingness_method, 
    leverage=PC.leverage,
    robdist=PC.robdist,
    robdist_subset=PC.robdist_subset)
  id_out.fun <- switch(outlyingness_method, 
    leverage=id_out.leverage,
    robdist=id_out.robdist,
    robdist_subset=id_out.robdist_subset)

  N_ <- ncol(X)
  T_ <- nrow(X)
  if(N_ < T_){
    warning("Data matrix has more rows than columns. Check that observations
      are in rows and variables are in columns.\n")
  }

  # Center and scale the data robustly.
  # Do it here instead of calling scale_med to save memory.
  if(verbose){ print("Centering and scaling the data matrix.") }
  X <- t(X)
  # Center.
  X <- X - c(rowMedians(X, na.rm=TRUE))
  # Scale.
  mad <- 1.4826 * rowMedians(abs(X), na.rm=TRUE)
  zero_mad <- mad < TOL
  if(any(zero_mad)){
    if(all(zero_mad)){
      stop("Error: All voxels are zero-variance. \n")
    } else {
      warning(paste0("Warning: ", sum(zero_mad),
        " zero-variance voxels (out of ", length(zero_mad),
        "). These will be set to zero for estimation of the covariance.\n"))
    }
    mad[zero_mad] <- 1
  }
  X <- X/c(mad)
  X[zero_mad,] <- 0
  X <- t(X)
  rm(mad, zero_mad)

  ##################
  # DATA PROJECTIONS
  ##################

  # Compute the PC scores (and directions, if leverage images or PCATF are desired).
  solve_directions <- (lev_img_level > 0) | ("PCATF" %in% projection_methods)
  if(verbose){
    print(paste0("Computing the",
                 ifelse(("PCA_var" %in% projection_methods) | ("PCA_kurt" %in% projection_methods),
                  ifelse("PCATF" %in% projection_methods, " normal and trend-filtered", ""),
                  ifelse("PCATF" %in% projection_methods, " trend-filtered", "INTERNAL ERROR")),
                 " PC scores",
                  ifelse(solve_directions, " and directions", ""), "."))
  }
  if(solve_directions){
    X.svd <- svd(X)
    rm(X)
  } else {
    # Avoid computing U matrix to conserve memory.
    XXt <- tcrossprod(X)
    rm(X)
    X.svd <- svd(XXt)
    rm(XXt)
    X.svd$d <- sqrt(X.svd$d)
    X.svd$v <- NULL
  }

  # Compute PCATF, if requested.
  if("PCATF" %in% projection_methods){
    X.svdtf <- do.call(
      PCATF, 
      c(list(X=X, X.svd=X.svd,
             K=choose_PCs.variance(X.svd, max_keep=NULL, min_keep=NULL),
             solve_directions=solve_directions),
      PCATF_kwargs))
    # The PC directions were needed to compute PCATF. If leverage images 
    #   are not wanted, we can now delete the directions to save space.
    if(lev_img_level > 0){ X.svd$v <- NULL }
  }
  gc()

  # Remove PCs with constant scores (often present in PCATF).
  remove_const_PCs <- function(the_svd, PC_name){
    zero_var <- apply(the_svd$u, 2, var) < TOL
    if(any(zero_var)){
      if(all(zero_var)){
        stop(paste("Error: All", PC_name, "scores are zero-variance.\n"))
      }
      warning(paste("Warning:", sum(zero_var), PC_name,
        "scores are zero-variance. Removing these PCs."))
      the_svd$u <- the_svd$u[,!zero_var]
      the_svd$d <- the_svd$d[!zero_var]
      the_svd$v <- the_svd$v[,!zero_var]
    }
  }
  X.svd <- remove_const_PCs(X.svd, 'PC')
  X.svdtf <- remove_const_PCs(X.svdtf, 'trend-filtered PC')

  # Choose which PCs to retain for each projection.
  projections <- setNames(vector("list", length(projection_methods)), projection_methods)
  for(i in 1:length(projections_methods)){
    projection_method = projections_methods[i]

    if(verbose){
      print(switch(projection_method,
        "PCA_var" = "Identifying the PCs with high varaince.",
        "PCA_kurt" = "Identifying the PCs with high kurtosis.",
        "PCATF" = "Identifying the trend-filtered PCs with high varaince."
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
      if(choose_PCs == "kurtosis"){
        choose_PCs_kwargs$kurt_quantile <- kurt_quantile
        choose_PCs_kwargs$detrend <- kurt_detrend
      }
      chosen_PCs <- do.call(choose_PCs.fun, choose_PCs_kwargs)
      chosen_PCs <- chosen_PCs[order(chosen_PCs)] # kurtosis order =/= index order
    }

    projections[projection_method] = list(
      indices = chosen_PCs,
      svd = ifelse(projection_method=="PCATF", X.svdtf, X.svd)
    )
    projections[projection_method]$svd$u <- projections[projection_method]$svd$u[,chosen_PCs]
    projections[projection_method]$svd$d <- projections[projection_method]$svd$d[chosen_PCs]
  }

  ##########################################
  # OUTLIER MEASUREMENTS AND LEVERAGE IMAGES
  ##########################################

  # Perform the rest of the steps for each method combination.
  outlier_meas <- outlier_cutoffs <- outlier_flags <- inMCD <- lev_imgs <- setNames(
    vector("list", length(method_combos)), method_combos)
  for(i in 1:length(method_combos)){
    if(verbose){
      print(paste0("Method ", method_combos, ":"))
    }
    method_combo = method_combos[i]
    method_split = unlist(strsplit(method_combo, "__"))
    projection_method = method_split[1]
    outlyingness_method = method_split[2]
    projection = projections[projection_method]

    # Adjust PC number if using robust distance (subset).
    if(outlyingness_method %in% c("robdist", "robdist_subset")){
      # Let q = N_/T_ (ncol(U)/nrow(U)). robustbase::covMcd()
      #  requires q <= approx. 1/2 for computation of the MCD covariance estimate.
      #  Higher q will use more components for estimation, thus retaining a
      #  higher resolution of information. Lower q will have higher breakdown
      #  points, thus being more resistant to outliers. (The maximal breakdown
      #  value is (N_ - T_ + 2)/2.) Here, we select q = 1/3 to yield a breakdown
      #  value of approx. 1/3. Since the subset method splits T_ into thirds, it
      #  must further reduce N_ by 1/3.
      q <- 1/3
      max_keep = ifelse(outlyingness_method=="robdist", T_*q, T_*q/3)
      if(max_keep > length(projection$indices)){
        projection$indices = projection$indices[1:max_keep]
        projection$svd$u = projection$svd$u[,1:max_keep]
        projection$svd$d = projection$svd$u[1:max_keep]
        if(solve_directions){
          projection$svd$v = projection$svd$v[,1:max_keep]
        }
      }
    }

    # Compute outlyingness measure.
    outlyingness_method.fun <- switch(
      "leverage" = id_out.leverage,
      "robdist" = id_out.robdist,
      "robdist_subset" = id_out.robdist_subset,
    )
    measure <- outlyingness_method.fun(projection$svd$u)
    if(outlyingness_method %in% c("robdist_subset", "robdist")){
      inMCD <- measure$inMCD
      Fparam <- measure$Fparam
      measure <- measure$robdist
    }
    outlier_meas[method_combo] <- measure
    
    if(outlyingness_method %in% c("robdist_subset", "robdist")){
      inMCD[method_combo] <- inMCD
    }
    
    # Identify outliers.
    if(id_out){
      if(verbose){ print("...identifying outliers.") }
      if(outlyingness_method == "leverage"){
        id_out_kwargs <- list(leverage=measure)
      }
      if(outlyingness_method %in% c("robdist_subset", "robdist")){
        id_out_kwargs <- list(distance=measure, inMCD=inMCD, Fparam=Fparam)
      }
      out <- do.call(id_out.fun, id_out_kwargs)
      outlier_cutoffs[method_combo] <- out$cutoffs
      outlier_flags[method_combo] <- out$outliers
    }

    # Make leverage images.
    if(lev_img_level > 0){
      if(sum(out$outliers[,lev_img_level]) > 0){
        if(verbose){
          print(paste0(
            "...outliers detected at level ", lev_img_level, " (",
            colnames(out$outliers)[lev_img_level], "). Computing leverage images."))
        }
        lev_imgs[method_combo] <- leverage_images(projection$svd, which(out$outliers[,lev_img_level]))
      } else {
        if(verbose){
          print(paste0(
            "...no leverage images: clever did not find any outliers at level ",
            lev_img_level, " (", colnames(out$outliers)[lev_img_level], ")."))
        }
        lev_imgs[method_combo] <- list(mean=NULL, top=NULL, top_dir=NULL)
      }
    }
  }

  #######################
  # DVARS & FORMAT OUTPUT
  #######################

  # Conserve memory.
  projections <- NULL
  gc()

  if(DVARS){
    DVARS <- compute_DVARS(X)
  } else {
    DVARS <- NULL
  }

  # Organize the output.
  if(verbose){ print("Done! Organizing results.") }
  result <- list(params=NULL, projections=projections, DVARS=DVARS,
                 outlier_meas=outlier_meas, outlier_cutoffs=outlier_cutoffs,
                 outlier_flags=outlier_flags, inMCD=inMCD, lev_imgs=lev_imgs)
  result$params <- list(
    projection_methods = projection_methods,
    outlyingness_methods = outlyingness_methods,
    DVARS = DVARS,
    PCATF_kwargs = PCATF_kwargs,
    kurt_quantile = kurt_quantile,
    kurt_detrend = kurt_detrend,
    id_out = id_out,
    lev_img_level = lev_img_level)

  class(result) <- c("clever", class(result))

  return(result)
}
