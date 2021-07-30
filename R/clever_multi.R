#' Compare multiple leverage measures with \code{clever_multi}
#' 
#' Calculates leverage to identify outliers in high-dimensional data. 
#'  Can compute multiple kinds of leverage at once. 
#' 
#' @inheritSection clever_order_of_operations Order of operations
#' @inheritParams clever_Params
#' @param projection Leverage works by projecting the data onto directions likely to 
#'  contain outlier information. Choose at least one of the following:
#' 
#'  \describe{
#'    \item{\code{"PCA"}}{PCA using the top \eqn{k} PCs.}
#'    \item{\code{"PCA_kurt"}}{PCA using the high-kurtosis PCs among the top \eqn{k}.}
#'    \item{\code{"PCA2"}}{PCA using the top \eqn{k2} PCs.}
#'    \item{\code{"PCA2_kurt"}}{PCA using the high-kurtosis PCs among the top \eqn{k2}.}
#'    \item{\code{"PCATF"}}{PCATF using the top \eqn{k} trend-filtered PCs.}
#'    \item{\code{"PCATF_kurt"}}{PCATF using the high-kurtosis trend-filtered PCs among the top \eqn{k}.}
#'    \item{\code{"PCATF2"}}{PCATF using the top \eqn{k2} trend-filtered PCs.}
#'    \item{\code{"PCATF2_kurt"}}{PCATF using the high-kurtosis trend-filtered PCs among the top \eqn{k2}.}
#'    \item{\code{"ICA"}}{ICA using the top \eqn{k} ICs.}
#'    \item{\code{"ICA_kurt"}}{ICA using the high-kurtosis ICs among the top \eqn{k}.}
#'    \item{\code{"ICA2"}}{ICA using the top \eqn{k2} ICs.}
#'    \item{\code{"ICA2_kurt"}}{ICA using the high-kurtosis ICs among the top \eqn{k2}.}
#'  }
#' 
#'  where \eqn{k} is the number of PCs selected by PESEL, and \eqn{k2} is the number
#'  of PCs with above-average variance.
#'  
#'  Use \code{"all"} to use all projection methods. Default: \code{"PCA_kurt"}.
#' @return A \code{"clever_multi"} object, i.e. a list with components
#' \describe{
#'  \item{measure}{A \eqn{T \times P} data.frame of numeric leverage values for each of the P projections in \code{projection}.}
#'  \item{outlier_cutoff}{A \eqn{1 \times P} data.frame of numeric outlier cutoff values for each projection (\code{cutoff} times the median leverage).}
#'  \item{outlier_flag}{A \eqn{T \times P} data.frame of logical values where \code{TRUE} indicates suspected outlier presence.}
#'  \item{mask}{
#'    A length \eqn{P} numeric vector corresponding to the data locations in \code{X}. Each value indicates whether the location was masked:
#'    \describe{
#'      \item{1}{The data location was not masked out.}
#'      \item{-1}{The data location was masked out, because it had at least one \code{NA} or \code{NaN} value.}
#'      \item{-2}{The data location was masked out, because it was constant.}
#'    }
#'  }
#'  \item{PCA}{
#'    This will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{get_dirs}.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}
#'      \item{nPCs_PESEL}{The number of PCs selected by PESEL.}
#'      \item{nPCs_avgvar}{The number of above-average variance PCs.}
#'    }
#'    where \code{Q} is the number of PCs selected by PESEL or of above-average variance (or the greater of the two if both were used). 
#'    If PCA was not used, all entries except \code{nPCs_PESEL} and/or \code{nPCs_avgvar} will not be included, depending on which
#'    method(s) was used to select the number of PCs.
#'  }
#'  \item{PCATF}{
#'    If PCATF was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{get_dirs}}
#'    }
#'  }
#'  \item{ICA}{
#'    If ICA was used, this will be a list with components:
#'    \describe{
#'      \item{S}{The \eqn{P \times Q} source signals matrix. Included only if \code{get_dirs}} 
#'      \item{M}{The \eqn{T \times Q} mixing matrix.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating mixing scores of high kurtosis.}
#'    }
#'  }
#' }
#'
#' @importFrom pesel pesel
#' @importFrom robustbase rowMedians
#' @importFrom stats mad qnorm var setNames
#'
#' @keywords internal
#' 
#' @examples
#' n_voxels = 1e4
#' n_timepoints = 100
#' X = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' clev = fMRIscrub:::clever_multi(X)
clever_multi = function(
  X, projection = "PCA_kurt", 
  nuisance=cbind(1, dct_bases(nrow(X), 4)),
  center=TRUE, scale=TRUE, comps_mean_dt=FALSE, comps_var_dt=FALSE,
  kurt_quantile=.99, PCATF_kwargs=NULL,
  get_dirs=FALSE, full_PCA=FALSE,
  get_outliers=TRUE, cutoff=4,
  verbose=FALSE){

  # ----------------------------------------------------------------------------
  # Check arguments. -----------------------------------------------------------
  # ----------------------------------------------------------------------------

  # `X` ------------------------------------------------------------------------
  if (verbose) { cat("Checking for missing, infinite, and constant data.\n") }
  X <- as.matrix2(X); class(X) <- "numeric"
  V0_ <- ncol(X)
  X_NA_mask <- apply(X, 2, function(x){any(x %in% c(NA, NaN, -Inf, Inf))})
  if (any(X_NA_mask)) {
    if (all(X_NA_mask)) { stop("All data columns have at least one NA, NaN, or infinte value. None of these are allowed.\n") }
    warning(
      "Removing ", sum(X_NA_mask),
      " columns with at least one NA, NaN, or infinite value (out of ", ncol(X), ").\n"
    )
    X <- X[,!X_NA_mask,drop=FALSE]
  }
  X_const_mask <- apply(X, 2, is_constant)
  if (any(X_const_mask)) {
    if (all(X_const_mask)) { stop("All data columns are constant.\n") }
    warning(
      "Removing ", sum(X_const_mask), 
      " constant data columns (out of ", ncol(X), ").\n"
    )
    X <- X[,!X_const_mask,drop=FALSE]
  }
  T_ <- nrow(X); V_ <- ncol(X)
  if (T_ > V_) {
    warning(
      "Data matrix has more rows than columns. Check that observations ",
      "are in rows and variables are in columns.\n"
    )
  }

  # Create output structure. ---------------------------------------------------
  out <- list(
    measure = list(),
    measure_info = list(),
    outlier_cutoff = list(),
    outlier_flag = list(),
    mask = rep(1, V0_),
    PCA = NULL,
    PCATF = NULL,
    ICA = NULL
  )

  out$mask[X_NA_mask] <- -1
  out$mask[!X_NA_mask][X_const_mask] <- -2

  # `projection`----------------------------------------------------------------
  valid_projection_PESEL <- c("PCA", "PCA_kurt", "PCATF", "PCATF_kurt", "ICA", "ICA_kurt")
  valid_projection_avgvar <- c("PCA2", "PCA2_kurt", "PCATF2", "PCATF2_kurt", "ICA2", "ICA2_kurt")
  valid_projection <- c(valid_projection_PESEL, valid_projection_avgvar)
  if ("all" %in% projection) {
    projection <- valid_projection
  } else {
    projection <- unique(match.arg(projection, valid_projection, several.ok=TRUE))
  }
  base_projection <- unique(gsub("2", "", gsub("_kurt", "", projection)))

  # `out$measure_info` ---------------------------------------------------------
  m_info <- data.frame(name = valid_projection[valid_projection %in% projection])
  m_info$type <- "Leverage"
  m_info$projection <- gsub("2", "", gsub("_kurt", "", m_info$name))
  m_info$PESEL <- !grepl("2", m_info$name)
  m_info$kurt <- grepl("kurt", m_info$name)
  m_info$comps_mean_dt <- comps_mean_dt
  m_info$comps_var_dt <- comps_var_dt
  out$measure_info <- m_info

  # `nuisance`------------------------------------------------------------------
  do_nuisance <- !(is.null(nuisance) || isFALSE(nuisance) || identical(nuisance, 0))
  if (do_nuisance) { 
    nuisance <- check_design_matrix(nuisance, T_)
    design_const_mask <- apply(nuisance, 2, is_constant)
    if (!any(design_const_mask)) {
      if (!any(abs(apply(X, 2, mean)) < 1e-8)) {
        warning("No intercept detected in `design`, yet the data are not centered.")
      }
    }
  }

  # other arguments ------------------------------------------------------------
  center <- as.logical(center); stopifnot(isTRUE(center) || isFALSE(center))
  scale <- as.logical(scale); stopifnot(isTRUE(scale) || isFALSE(scale))
  if (isTRUE(comps_mean_dt)) {
    comps_mean_dt <- 4
  } else if (isFALSE(comps_mean_dt)) { 
    comps_mean_dt <- 0
  } else {
    comps_mean_dt <- as.numeric(comps_mean_dt)
    stopifnot(comps_mean_dt >= 0)
  }
  if (isTRUE(comps_var_dt)) {
    comps_var_dt <- 4
  } else if (isFALSE(comps_var_dt)) { 
    comps_var_dt <- 0
  } else {
    comps_var_dt <- as.numeric(comps_var_dt)
    stopifnot(comps_var_dt >= 0)
  }
  comps_dt <- (comps_mean_dt > 0) || (comps_var_dt > 0)
  kurt_quantile <- as.numeric(kurt_quantile)
  stopifnot(kurt_quantile >= 0 && kurt_quantile <= 1)
  if(!identical(PCATF_kwargs, NULL)){
    names(PCATF_kwargs) <- match.arg(
      names(PCATF_kwargs), c("lambda", "niter_max", "TOL", "verbose"),
      several.ok=TRUE)
    if(length(names(PCATF_kwargs)) != length(unique(names(PCATF_kwargs)))){
      stop("Duplicate PCATF_kwargs were given.")
    }
  }
  get_dirs <- as.logical(get_dirs); stopifnot(isTRUE(get_dirs) || isFALSE(get_dirs))
  full_PCA <- as.logical(full_PCA); stopifnot(isTRUE(full_PCA) || isFALSE(full_PCA))
  get_outliers <- as.logical(get_outliers); stopifnot(isTRUE(get_outliers) || isFALSE(get_outliers))
  cutoff <- as.numeric(cutoff); stopifnot(cutoff >= 0)
  verbose <- as.logical(verbose); stopifnot(isTRUE(verbose) || isFALSE(verbose))

  # ----------------------------------------------------------------------------
  # Nuisance regression followed by centering & scaling. -----------------------
  # ----------------------------------------------------------------------------

  if (do_nuisance) { 
    if (verbose) { cat("Performing nuisance regression.\n") }
    X <- nuisance_regression(X, nuisance)
  }

  if (center || scale) {
    if (verbose) { 
      action <- c("Centering", "Scaling", "Centering and scaling")[1*center + 2*scale]
      cat(action, "the data columns.\n")
    }

    # Center & scale here, rather than calling `scale_med`, to save memory.
    X <- t(X)
    if (center) { X <- X - c(rowMedians(X)) }
    if (scale) { X <- X / (1.4826 * rowMedians(abs(X))) }
    X <- t(X)
  }

  # ----------------------------------------------------------------------------
  # Make projection. -----------------------------------------------------------
  # ----------------------------------------------------------------------------

  # Compute PCA. (Even if only ICA is being used, since we need `nComps`)
  if (verbose) {
    cat(paste0(
      "Computing PCA.\n"
    ))
  }
  if (get_dirs || "PCATF" %in% base_projection) {
    out$PCA <- tryCatch(
      {
        svd(X)[c("u", "d", "v")]
      },
      error = function(cond) {
        message(cond)
        cat(
          "Trying `corpcor::fast.svd`. An error will occur if this package ",
          "is not available, in which case the package should be installed ",
          "and `clever` should be run again.\n"
        )
        if (!requireNamespace("corpcor", quietly = TRUE)) {
          stop("Package \"corpcor\" needed since `svd` failed. Please install it.", call. = FALSE)
        }
        return(corpcor::fast.svd(X)[c("u", "d", "v")])
      }
    )
    names(out$PCA) <- toupper(names(out$PCA))
  } else {
    # Conserve memory by using `XXt`.
    out$PCA <- tryCatch(
      {
        svd(tcrossprod(X))[c("u", "d", "v")]
      },
      error = function(cond) {
        message(cond)
        cat(
          "Trying `corpcor::fast.svd`. An error will occur if this package",
          "is not available, in which case the package should be installed",
          "and `clever` should be run again.\n"
        )
        if (!requireNamespace("corpcor", quietly = TRUE)) {
          stop("Package \"corpcor\" needed since `svd` failed. Please install it.", call. = FALSE)
        }
        return(corpcor::fast.svd(tcrossprod(X))[c("u", "d", "v")])
      }
    )
    names(out$PCA) <- toupper(names(out$PCA))
    out$PCA$D <- sqrt(out$PCA$D)
    out$PCA$V <- NULL
  }
  # Keep only the above-average variance/PESEL PCs (whichever is greater).
  maxK_PCA <- 1
  # [TO DO]: move outside if(PCA) so PCA isn't required for ICA+PESEL
  if (any(valid_projection_PESEL %in% projection)) {
    out$PCA$nPCs_PESEL <- with(set.seed(0), pesel::pesel(t(X), npc.max=ceiling(T_/2), method="homogenous")$nPCs)
    maxK_PCA <- max(maxK_PCA, out$PCA$nPCs_PESEL)
  }
  if (any(valid_projection_avgvar %in% projection)) {
    out$PCA$nPCs_avgvar <- max(1, sum(out$PCA$D^2 > mean(out$PCA$D^2)))
    maxK_PCA <- max(maxK_PCA, out$PCA$nPCs_avgvar)
  }
  # Identify which PCs have high kurtosis.
  if (any(c("PCA_kurt", "PCA2_kurt") %in% projection)) {
    out$PCA$highkurt <- high_kurtosis(out$PCA$U[, seq(maxK_PCA), drop=FALSE], kurt_quantile=kurt_quantile)
  }
  # [TO DO]: Resolve case where no PC has high kurtosis

  # Compute PCATF.
  if ("PCATF" %in% base_projection) {
    maxK_PCATF <- max(as.numeric(list(
      PCATF = out$PCA$nPCs_PESEL,
      PCATF_kurt = out$PCA$nPCs_PESEL,
      PCATF2 = out$PCA$nPCs_avgvar,
      PCATF2_kurt = out$PCA$nPCs_avgvar
    )[projection[grepl("PCATF", projection)]]))
    if (verbose) { cat("Computing PCATF.\n") }
    out$PCATF <- do.call(
      PCATF, 
      c(
        list(
          X=X, X.svd=out$PCA[c("U", "D", "V")], 
          K=maxK_PCATF, solve_directions=get_dirs
        ), 
        PCATF_kwargs
      )
    )
    out$PCATF$PC_exec_times <- NULL; out$PCATF$nItes <- NULL
    # V matrix from PCA no longer needed.
    if(!get_dirs){ out$PCA$V <- NULL }

    tf_const_mask <- apply(out$PCATF$u, 2, is_constant)
    if(any(tf_const_mask)){
      warning(
        "Warning: ", sum(tf_const_mask), " out of ", length(tf_const_mask),
        "trend-filtered PC scores are zero-variance.\n"
      )
    }
    names(out$PCATF)[names(out$PCATF) %in% c("u", "d", "v")] <- toupper(names(out$PCATF)[names(out$PCATF) %in% c("u", "d", "v")])
  }
  # Identify which trend-filtered PCs have high kurtosis.
  if (any(c("PCATF_kurt", "PCATF2_kurt") %in% projection)) {
    out$PCATF$highkurt <- high_kurtosis(out$PCATF$U, kurt_quantile=kurt_quantile)
  }

  # Compute ICA
  if (any(c("ICA", "ICA2") %in% base_projection)) {
    maxK_ICA <- max(as.numeric(list(
      ICA = out$PCA$nPCs_PESEL,
      ICA_kurt = out$PCA$nPCs_PESEL,
      ICA2 = out$PCA$nPCs_avgvar,
      ICA2_kurt = out$PCA$nPCs_avgvar
    )[projection[grepl("ICA", projection)]]))
    if (verbose) { cat("Computing ICA.\n" ) }
    if (!requireNamespace("ica", quietly = TRUE)) {
      stop("Package \"ica\" needed to compute the ICA. Please install it.", call. = FALSE)
    }
    out$ICA <- with(set.seed(0), ica::icaimax(t(X), maxK_ICA, center=FALSE))[c("S", "M")]
    # Issue due to rank.
    if (ncol(out$ICA$M) != maxK_ICA) {
      cat("Rank issue with ICA: adding constant zero columns.\n")
      K_missing <- maxK_ICA - ncol(out$ICA$M)
      out$ICA$M <- cbind(out$ICA$M, matrix(0, nrow=nrow(out$ICA$M), ncol=K_missing))
      out$ICA$S <- cbind(out$ICA$S, matrix(0, nrow=nrow(out$ICA$S), ncol=K_missing))
    }

    if (any(c("ICA_kurt", "ICA2_kurt") %in% projection)) {
      out$ICA$highkurt <- high_kurtosis(out$ICA$M, kurt_quantile=kurt_quantile)
    }

    if(!get_dirs){ out$ICA$S <- NULL }
  }

  # Remove PCA information if only ICA is being used.
  # Do this after PCA info was given to PCATF
  if (!("PCA" %in% base_projection)) {
    out$PCA$U <- out$PCA$D <- out$PCA$V <- NULL
  } else {
    # Remove smaller PCs.
    if (!full_PCA) {
      out$PCA$U <- out$PCA$U[, seq(maxK_PCA), drop=FALSE]
      out$PCA$D <- out$PCA$D[seq(maxK_PCA), drop=FALSE]
      if (!is.null(out$PCA$V)) { 
        out$PCA$V <- out$PCA$V[, seq(maxK_PCA), drop=FALSE]
      }
    }
  }

  if (comps_dt) {
    if (!is.null(out$PCA$U)) { 
      out$PCA$U_dt <- apply(out$PCA$U, 2, rob_scale, center=comps_mean_dt, scale=comps_var_dt)
      out$PCA$highkurt_dt <- high_kurtosis(out$PCA$U_dt, kurt_quantile=kurt_quantile)
    }
    if (!is.null(out$PCATF$U)) { 
      out$PCATF$U_dt <- apply(out$PCATF$U, 2, rob_scale, center=comps_mean_dt, scale=comps_var_dt)
      out$PCATF$highkurt_dt <- high_kurtosis(out$PCATF$U_dt, kurt_quantile=kurt_quantile)
    }
    if (!is.null(out$ICA$M)) { 
      out$ICA$M_dt <- apply(out$ICA$M, 2, rob_scale, center=comps_mean_dt, scale=comps_var_dt)
      out$ICA$highkurt_dt <- high_kurtosis(out$ICA$M_dt, kurt_quantile=kurt_quantile)
    }
  }

  rm(X); gc()

  # ----------------------------------------------------------------------------
  # Compute leverage. ----------------------------------------------------------
  # ----------------------------------------------------------------------------

  if (verbose) { cat("Computing leverage.\n") }

  for (ii in seq(length(projection))) {
    proj_ii <- projection[ii]
    base_ii <- gsub("2", "", gsub("_kurt", "", proj_ii))
    scores_ii <- paste0(
      ifelse(grepl("ICA", proj_ii), "M", "U"), 
      ifelse(!isFALSE(comps_dt), "_dt", "")
    )
    highkurt_ii <- paste0("highkurt", ifelse(comps_dt, "_dt", ""))

    # Make projection.
    Comps_ii <- switch(proj_ii,
      PCA = seq(out$PCA$nPCs_PESEL),
      PCA_kurt = which(out$PCA[[highkurt_ii]][seq(out$PCA$nPCs_PESEL)]),
      PCA2 = seq(out$PCA$nPCs_avgvar),
      PCA2_kurt = which(out$PCA[[highkurt_ii]][seq(out$PCA$nPCs_avgvar)]),
      PCATF = seq(out$PCA$nPCs_PESEL),
      PCATF_kurt = which(out$PCATF[[highkurt_ii]][seq(out$PCA$nPCs_PESEL)]),
      PCATF2 = seq(out$PCA$nPCs_avgvar),
      PCATF2_kurt = which(out$PCATF[[highkurt_ii]][seq(out$PCA$nPCs_avgvar)]),
      ICA = seq(out$PCA$nPCs_PESEL),
      ICA_kurt = which(out$ICA[[highkurt_ii]][seq(out$PCA$nPCs_PESEL)]),
      ICA2 = seq(out$PCA$nPCs_avgvar),
      ICA2_kurt = which(out$ICA[[highkurt_ii]][seq(out$PCA$nPCs_avgvar)])
    )

    if (grepl("kurt", proj_ii) && length(Comps_ii) < 1) {
      result_ii <- list(
        meas = rep(0, T_), 
        cut = NA, 
        flag = rep(FALSE, T_)
      )
    } else {
      Comps_ii <- out[[base_ii]][[scores_ii]][, Comps_ii, drop=FALSE]
      # Compute leverage.
      result_ii <- out_measures.leverage(Comps=Comps_ii, median_cutoff=cutoff)
    }
    out$measure[[proj_ii]] <- result_ii$meas
    if (get_outliers) {
      out$outlier_cutoff[[proj_ii]] <- result_ii$cut
      out$outlier_flag[[proj_ii]] <- result_ii$flag
    }
  }

  # ----------------------------------------------------------------------------
  # Format output. -------------------------------------------------------------
  # ----------------------------------------------------------------------------

  out$measure <- as.data.frame(out$measure)
  if (length(out$outlier_cutoff) > 0) {
    out$outlier_cutoff <- as.data.frame(out$outlier_cutoff)
  } else {
    out$outlier_cutoff <- NULL
  }
  if (length(out$outlier_flag) > 0) {
    out$outlier_flag <- as.data.frame(out$outlier_flag)
  } else {
    out$outlier_flag <- NULL
  }

  out <- out[!vapply(out, is.null, FALSE)]
 
  structure(out, class="clever_multiLev")
}