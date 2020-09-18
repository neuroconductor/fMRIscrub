#' Computes PCA leverage.
#'
#' Computes the leverage of each observation, in the PC score matrix (U).
#'  Optionally can identify the outliers.
#'
#' @param U An n x Q matrix of PC scores.
#' @param median_cutoff The outlier cutoff, in multiples of the median leverage.
#'  Default: \code{NA} (do not compute outliers).
#' 
#' @return A list with entries \code{"meas"} (the leverage values), 
#'  \code{"cut"} (the leverage cutoff value) and 
#'  \code{"flag"} (logical vector indicating the outliers). If 
#'  \code{!is.null(median_cutoff)}, all entries except \code{"meas"} are omitted.
#'  
#' @export
out_measures.leverage <- function(U, median_cutoff=NA){
  # Below line: same as diag(U %*% t(U)), but faster.
  lev <- apply(U, 1, function(x){sum(x*x)})
  out <- list(meas=lev)
  if (!is.null(median_cutoff)){
    out$cut <- median_cutoff * median(lev)
    out$flag <- out$meas > out$cut
  }
  out
}

#' Computes MCD distances.
#'
#' Computes robust minimum covariance determinant (MCD) distances across
#'  the observations (rows). The MCD method selects a subset of h observations
#'  whose covariance matrix has minimum determinant across all subsets of size
#'  h. The MCD distances are Mahalanobis distances using the estimates of
#'  center (mean) and scale (covariance matrix) based on that subset.
#'
#' @param U An n x Q matrix of PC scores.
#' @param quantile_cutoff The F-distribution quantile cutoff. Default: 
#'  \code{NA} (do not compute outliers).
#' 
#' @return A list with entries
#' \describe{
#'   \item{"meas"}{A vector of length n of with the robust distance estimate
#'    for each observation.}
#'   \item{"info"}{A list with entries "inMCD", "outMCD_scale", and "Fparam"}
#'   \item{"cut"}{The robust distance cutoff value}
#'   \item{"flag"}{Logical vector indicating the outliers}
#' }
#' 
#' If \code{is.null(quantile_cutoff)} the latter two elements are omitted.
#'
#' @importFrom MASS cov.mcd
#'
#' @export
out_measures.robdist <- function(U, quantile_cutoff=NA){ 
  n <- nrow(U)
  Q <- ncol(U)
  
  best <- c(cov.mcd(U)$best)
  inMCD <- 1:n %in% best
  U_in <- matrix(U[best,], ncol=Q) # observations that are used for the MCD estimates calculation
  xbar_star <- colMeans(U_in) # MCD estimate of mean
  U_ins <- scale(U_in, center = TRUE, scale = FALSE) 
  nU <- nrow(U_ins)
  S_star <- (t(U_ins) %*% U_ins)/(nU-1) # MCD estimate of covariance
  RD <- (apply(U, 1, function(x) t(x-xbar_star) %*% solve(S_star) %*% (x-xbar_star)))
  # Scale left-out observations to follow F-distribution.
  Fparam <- fit.F(Q, n, sum(inMCD))
  Fparam <- c(Fparam$c, Fparam$m, Fparam$df[1], Fparam$df[2])
  names(Fparam) <- c("c", "m", "df1", "df2")
  outMCD_scale <- Fparam["c"] * (Fparam["m"] - Q + 1) / (Q * Fparam["m"]) 
  
  out <- list(
    meas=RD, 
    info = list(inMCD=inMCD, outMCD_scale=outMCD_scale, Fparam=Fparam)
  )

  if (!is.null(quantile_cutoff)) {
    out$cut <- qf(p=quantile_cutoff, df1=Fparam["df1"], df2=Fparam["df2"])
    out$flag <- ifelse(inMCD, FALSE, RD * outMCD_scale > out$cut)
  }

  out
}

#' Identify outliers based on bootstrap robust distance.
#' 
#' @param U N x P Dependent data
#' @param R_true N x N correlation matrix. If \code{NULL} (default), the identity
#'  matrix wil lbe used
#' @param boot_sample Number of bootstrap samples to compute
#' @param quantile_cutoff The F-distribution quantile cutoff. Default: 
#'  \code{NA} (do not compute outliers).
#' 
#' @return List of "meas" (original RD values), "cut" (P cutoff values, one for 
#'  each observation) and "flag" (P indicators of outlyingness). If \code{is.null(quantile_cutoff)},
#'  the latter two entries are omitted.
#' 
#' @importFrom expm sqrtm
#' 
out_measures.robdist_bootstrap <- function(U, R_true=NULL, boot_sample=1000, quantile_cutoff=NA){
  if (is.null(R_true)) { R_true <- diag(nrow(U)) }

  x < - induc_indep(U, R_true)
  Y_tilde <- x$indep_data
  R_sqrt <- x$R_sqrt
  n <- nrow(Y_tilde)
  p <- ncol(Y_tilde)

  out1 <- out_measures.robdist(Y_tilde, quantile_cutoff=quantile_cutoff)
  best <- which(out1$inMCD)
  h <- length(best)
  
  # Obtain "notout"
  notbest <- setdiff(1:n, best)
  id_out <- which(out1$flag)
  notout <- setdiff(notbest, id_out)
  
  Y_tilde_in <- Y_tilde[best,]
  xbar_star <- colMeans(Y_tilde_in) # MCD estimate of mean , length of p
  Y_tilde_ins <- scale(Y_tilde_in, center = TRUE, scale = FALSE )
  nY_tilde <- nrow(Y_tilde_ins)
  S_star <- (t(Y_tilde_ins) %*% Y_tilde_ins)/(nY_tilde-1) # MCD estimate of covariance
  invcov <- solve(S_star) # dim p by p
  
  # Robust Distance of the Original Data
  RD_orig <- apply(U, 1, function(x) t(x-xbar_star) %*% invcov %*% (x-xbar_star))
  out <- list(meas=RD_orig)

  if (!is.null(quantile_cutoff)){
    # pre-compute these to use for quickly computing robust distances for many observations simultaneously with matrix operations
    invcov_sqrt <- sqrtm(invcov)
    xbar_star_mat <- matrix(xbar_star, nrow=n, ncol=p , byrow = TRUE)
    
    # QUANTILE ADJUSTMENT for the "mostly notbest" observations
    gamma <- length(notbest)/ n
    alpha <- 1- quantile_cutoff
    adj_alpha <- alpha * gamma
    ### BOOTSTRAP STEP
    Y_tilde_boot <- array(NA, dim=c(n,p))
    RD_boot <- matrix(NA, n, boot_sample)
    for(b in 1:boot_sample){ # start loop over bootstrap samples
      Y_tilde_boot <- 0 * Y_tilde_boot
      best_boot <- sample(best, h, replace = T)
      notout_boot <- sample(notout, (n-h), replace=T)
      Y_tilde_boot[best,] <- Y_tilde[best_boot,]
      Y_tilde_boot[notbest,] <- Y_tilde[notout_boot,]
      Y_boot_b <- R_sqrt %*% Y_tilde_boot #RE-INDUCING DEPENDENCE
      #RD_boot[,b] <- apply(Y_boot_b,1, function(x) t(x-xbar_star) %*% invcov %*% (x-xbar_star))
      #following two lines are alternative way to get RDs
      temp <- (Y_boot_b - xbar_star_mat) %*% invcov_sqrt
      RD_boot[,b] <- rowSums(temp * temp)
    } # end loop over bootstrap samples
    RD_scaled <- out1$outMCD_scale * RD_boot # scale is from Hardin & Rocke's approach (2005)
    
    out$cut <- t(apply(RD_scaled[notbest,], 1, quantile, probs=c(1-adj_alpha)))
    out$flag <- rep(FALSE, n)
    out$flag[notbest] <- (RD_orig[notbest] > out$cut)
  }

  out
}
