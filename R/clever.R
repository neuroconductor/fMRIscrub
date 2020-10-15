#' Identify outliers with \code{clever}
#' 
#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' \code{clever} will use all combinations of the requested projection and 
#'  out_meas methods that make sense. For example, if  
#'  \code{projection=c("PCATF", "PCA_var", "PCA_kurt")} and 
#'  \code{out_meas=c("leverage", "robdist")} then these five
#'  combinations will be used: PCATF with leverage, PCA + variance with 
#'  leverage, PCA + variance with robust distance, PCA + kurtosis with leverage,
#'  and PCA + kurtosis with robust distance. Each method combination will yield 
#'  its own out_meas time series.
#' 
#' @param X Numerical data matrix. Should be wide (N observations x P variables, 
#'  \eqn{N >> P}).
#' @param projection Character vector indicating the projection methods
#'  to use. Choose at least one of the following: \code{"PCA_var"} for 
#'  PCA + variance, \code{"PCA_kurt"} for PCA + kurtosis, and \code{"PCATF"} for
#'  PCA Trend Filtering + variance. Or, use \code{"all"} to use all projection 
#'  methods. Default: \code{c("PCA_kurt")}.
#' @param out_meas Character vector indicating the outlyingness measures to 
#'  compute. Choose at least one of the following: \code{"leverage"} for 
#'  leverage, \code{"robdist"} for robust distance, or \code{"robdist_bootstrap"}
#'  for robust distane bootstrap (not implemented yet). Or, use \code{"all"} 
#'  to use all methods. Default: \code{c("leverage")}.
#' @param DVARS Should DVARS (Afyouni and Nichols, 2017) be computed too? Default 
#'  is \code{TRUE}.
#' @param detrend Detrend the PCs before measuring kurtosis or computing 
#'  leverage or robust distance? Default: \code{TRUE}.
#' 
#'  Detrending is highly recommended for time-series data, especially if there 
#'  are many time points or evolving circumstances affecting the data. There are
#'  two reasons: first, temporal trends induce positive or negative kurtosis, 
#'  contaminating the connection between high kurtosis and outlier presence. 
#'  Second, trends tend to reduce the size of the in-MCD subset for the robust 
#'  distance method causing many false positives.
#'  
#'  Detrending should not be used with non-time-series data because the 
#'  observations are not temporally related.
#' 
#'  In addition to \code{TRUE} and \code{FALSE}, a third option \code{kurtosis}
#'  can be used to only detrend the PCs for the purpose of measuring kurtosis, 
#'  and not for the actual outlyingness measurement.
#' 
#'  This option will not affect the PCATF projection, which is never detrended.
#' @param PCATF_kwargs Named list of arguments for PCATF projection method.
#'  Only applies if \code{("PCATF" \%in\% projection)}.
#' 
#'  Valid entries are: 
#'  
#'  \describe{
#'    \item{K}{maximum number of PCs to compute (Default: \code{1000})}
#'    \item{lambda}{trend filtering parameter (Default: \code{0.05})}
#'    \item{niter_max}{maximum number of iterations (Default: \code{1000})}
#'    \item{verbose}{Print updates? (Default: \code{FALSE})}
#'  }
#' @param kurt_quantile What cutoff quantile for kurtosis should be used? Only 
#'  applies if \code{("PCA_kurt" \%in\% projection)}. Default: \code{0.95}.
#' @param id_outliers Should the outliers be identified? Default: \code{TRUE}.
#' @param lev_cutoff The outlier cutoff value for leverage, as a multiple of the
#'  median leverage. Only used if 
#'  \code{"leverage" \%in\% projection} and \code{id_outliers}. Default: 
#'  \code{4}, or \eqn{4 * median}.
#' @param rbd_cutoff The outlier cutoff quantile for MCD distance. Only used if 
#'  \code{"robdist" \%in\% projection} and \code{id_outliers}. Default: 
#'  \code{0.9999}, for the \eqn{0.9999} quantile.
#' 
#'  The quantile is computed from the estimated F distribution.
#' @param R_true The N x N correlation matrix, if known. Used for the bootstrap
#'  robust distance method.
#' @param lev_images Should leverage images be computed? If \code{FALSE} memory
#'  is conserved. Default: \code{FALSE}.
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#'  \item{params}{A list of all the arguments used.}
#'  \item{projections}{
#'    \describe{
#'      \item{PC_var}{
#'        \describe{
#'          \item{indices}{The indices retained from the original SVD 
#'            projection to make the variance-based PC projection.} 
#'          \item{PCs}{The PC projection.}  
#'        }
#'      }
#'      \item{PC_kurt}{
#'        \describe{
#'          \item{indices}{The indices retained from the original SVD 
#'            projection to make the kurtosis-based PC projection. They are 
#'            ordered from highest kurtosis to lowest kurtosis.}  
#'          \item{PCs}{The PC projection. PCs are ordered in the standard
#'            way, from highest variance to lowest variance, instead of by 
#'            kurtosis.}  
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
#'      \item{PC_var__leverage}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{PC_kurt__leverage}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{PCATF__leverage}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{PC_var__robdist}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{PC_kurt__robdist}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{DVARS_DPD}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'      \item{DVARS_ZD}{Logical vector idnicating whether each observation surpasses the outlier cutoff.}
#'    }
#'  }
#'  \item{robdist_info}{
#'    \describe{
#'      \item{PC_var__robdist}{
#'        \describe{
#'          \item{inMCD}{Logical vector indicating whether each observation was in the MCD estimate.}
#'          \item{outMCD_scale}{The scale for out-of-MCD observations.}
#'          \item{Fparam}{Named numeric vector: c, m, df1, and df2.}
#'        }
#'      }
#'      \item{PC_var__robdist}{
#'        \describe{
#'          \item{inMCD}{Logical vector indicating whether each observation was in the MCD estimate.}
#'          \item{outMCD_scale}{The scale for out-of-MCD observations.}
#'          \item{Fparam}{Named numeric vector: c, m, df1, and df2.}
#'        }
#'      }
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
  projection = "PCA_kurt",
  out_meas = "leverage",
  DVARS = TRUE,
  detrend = TRUE,
  PCATF_kwargs = NULL,
  kurt_quantile = .95,
  id_outliers = TRUE,
  lev_cutoff = 4,
  rbd_cutoff = 0.9999,
  R_true = NULL,
  lev_images = FALSE,
  verbose = FALSE) {

  # ----------------------------------------------------------------------------
  # Check arguments. -----------------------------------------------------------
  # ----------------------------------------------------------------------------

  # Define the valid `projection` and `out_meas`, and their valid combos.
  all_projection <- c("PCA_var", "PCA_kurt", "PCATF")
  all_out_meas <- c("leverage", "robdist")#, "robdist_bootstrap")
  all_valid_methods <- c(
    "PCA_var__leverage", "PCA_kurt__leverage", 
    "PCA_var__robdist", "PCA_kurt__robdist", 
    #"PCA_var__robdist_bootstrap", "PCA_kurt__robdist_bootstrap",
    "PCATF__leverage"
  )

  # Define the cutoff value for detecting zero variance/MAD voxels
  TOL <- 1e-8

  # Check arguments.
  if(!is.matrix(X)){ X <- as.matrix(X) }
  class(X) <- "numeric"

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
  stopifnot(is.TRUEorFALSE(DVARS))
  stopifnot(is.TRUEorFALSE(id_outliers))
  stopifnot(is.TRUEorFALSE(lev_images))
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
  
  if((lev_images) && (!id_outliers)){
    stop(
      "Invalid argument: computing leverage images requires\
      `id_outliers==TRUE`."
    )
  }

  stopifnot(is.numeric(lev_cutoff)); stopifnot(lev_cutoff > 0)

  stopifnot((rbd_cutoff > 0) & (rbd_cutoff < 1))

  # Get data dimensions.
  Npre_ <- ncol(X)
  T_ <- nrow(X)
  if(Npre_ < T_){
    warning(
      "Data matrix has more rows than columns. Check that observations\
      are in rows and variables are in columns."
    )
  }

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
  if(id_outliers){
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

    if(id_outliers){
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
      if (id_outliers) { cat(":") }
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
    if (id_outliers) {
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
    if (id_outliers) {
      outlier_cutoffs[[method_ii]] <- out_ii$cut
      outlier_flags[[method_ii]] <- out_ii$flag
    }
    if (out_ii_name == "robdist") {
      robdist_info[[method_ii]] <- out_ii$info
    }
      
    # --------------------------------------------------------------------------
    # Make leverage images.-----------------------------------------------------
    # --------------------------------------------------------------------------

    if (id_outliers) {
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
  if (id_outliers) {
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
    id_outliers = id_outliers,
    lev_cutoff = lev_cutoff,
    rbd_cutoff = rbd_cutoff,
    lev_images = lev_images,
    verbose = verbose
  )

  structure(result, class="clever")
}