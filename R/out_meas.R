#' Computes PCA leverage.
#'
#' Computes the leverage of each observation, based on its PC scores.
#'
#' @param U An n x Q matrix of PC scores.
#'
#' @return A vector of length n with the leverage of each observation.
#' @export
PC.leverage <- function(U){
  #diag(U %*% t(U))
  # Below line yields identical results as above line, but it's faster.
  apply(U, 1, function(x){sum(x*x)})
}

#' Computes MCD distances.
#'
#' Computes robust minimum covariance determinant (MCD) distances across
#'  the observations (rows).  The MCD method selects a subset of h observations
#'  whose covariance matrix has minimum determinant across all subsets of size
#'  h. The MCD distances are Mahalanobis distances using the estimates of
#'  center (mean) and scale (covariance matrix) based on that subset.
#'
#'
#' @param U An n x Q matrix of PC scores.
#'
#' @return A list with components
#' \describe{
#'   \item{robdist}{A vector of length n of with the robust distance estimate
#'    for each observation.}
#'  \item{inMCD}{A vector of length n indicating if each observation is within
#'    the MCD subset.}
#'  \item{Fparam}{The estimated parameters of the F distribution of MCD
#'    distances.}
#' }
#'
#' @importFrom MASS cov.mcd
#'
#' @export
PC.robdist <- function(U){ 
  n <- nrow(U)
  Q <- ncol(U)
  
  best <- c(cov.mcd(U)$best)
  inMCD <- 1:n %in% best
  U_in <- matrix(U[best,], ncol=Q) # observations that are used for the MCD estimates calculation
  xbar_star <- colMeans(U_in) # MCD estimate of mean
  U_ins <- scale(U_in, center = TRUE, scale = FALSE) # Damon: replace with U_ins - xbar_star?
  nU <- nrow(U_ins)
  S_star <- (t(U_ins) %*% U_ins)/(nU-1) # MCD estimate of covariance
  mah <- (apply(U, 1, function(x) t(x-xbar_star) %*% solve(S_star) %*% (x-xbar_star)))
  # Scale left-out observations to follow F-distribution.
  Fparam <- fit.F(Q, n, sum(inMCD))
  outMCD_scale <- Fparam$c * (Fparam$m - Q + 1) / (Q * Fparam$m) 
  #mah[!inMCD] <- outMCD_scale*mah[!inMCD] # Now, return scale and mah separately.
  
  result <- list(mah, inMCD, outMCD_scale, Fparam)
  names(result) <- c("mah", "inMCD", "outMCD_scale", "Fparam")
  return(result)
}

induc_indep <- function(dep_data, corr_matrix){
  U <- svd(R_true)$u
  V <- svd(R_true)$v
  D_invsqrt <- diag(1/sqrt(svd(R_true)$d))
  R_invsqrt <- U %*% D_invsqrt %*% t(V) # inverse square root of correlation matrix
  Y_tilde <- R_invsqrt %*% dep_data # n by p
  
  # for re-inducing dependence later
  D_sqrt <- diag(sqrt(svd(R_true)$d))
  R_sqrt <- U %*% D_sqrt %*% t(V)
  
  result <- list(Y_tilde,R_sqrt)
  names(result) <- c('indep_data','R_sqrt')
  return(result)
}

id_out_boot_robdist <- function(dep_data,boot_sample){
  Y_tilde <- induc_indep(dep_data,R_true)$indep_data
  R_sqrt <- induc_indep(dep_data,R_true)$R_sqrt
  scale <- PC.robdist(Y_tilde)$scale
  n <- dim(Y_tilde)[1]
  p <- dim(Y_tilde)[2]
  B <- boot_sample
  best <- cov.mcd(Y_tilde, nsamp="best")$best # AKA inMCD
  h <- length(best)
  
  # start to obtain "notout"
  notbest<- setdiff(1:n, best)
  out <- id_out.robdist(PC.robdist(Y_tilde)$robdist, PC.robdist(Y_tilde)$inMCD,PC.robdist(Y_tilde)$Fparam)
  id_out <- c(which(out$outliers[,3] ==T))
  notout <- setdiff(notbest, id_out)
  # end to obtain "notout"
  
  Y_tilde_in <- Y_tilde[best,]
  xbar_star <- colMeans(Y_tilde_in) # MCD estimate of mean , length of p
  Y_tilde_ins <- scale(Y_tilde_in, center = TRUE, scale = FALSE )
  nY_tilde <- nrow(Y_tilde_ins)
  S_star <- (t(Y_tilde_ins) %*% Y_tilde_ins)/(nY_tilde-1) # MCD estimate of covariance
  invcov <- solve(S_star) # dim p by p
  
  # pre-compute these to use for quickly computing robust distances for many observations simultaneously with matrix operations
  invcov_sqrt <- sqrtm(invcov)
  xbar_star_mat <- matrix(xbar_star, nrow=n, ncol=p , byrow = TRUE)
  
  # QUANTILE ADJUSTMENT for the "mostly notbest" observations
  gamma <- length(notbest)/ n
  alpha <- c(0.01, 0.001, 0.0001)
  adj_alpha <- alpha * gamma
  ### BOOTSTRAP STEP
  Y_tilde_boot <- array(NA,dim=c(n,p))
  RD_boot <- matrix(NA,n,B)
  for(b in 1:B){ # start loop over bootstrap samples
    Y_tilde_boot <- 0 * Y_tilde_boot
    best_boot <- sample(best, h, replace = T)
    notout_boot <- sample(notout, (n-h), replace=T)
    Y_tilde_boot[best,] <- Y_tilde[best_boot,]
    Y_tilde_boot[notbest,] <- Y_tilde[notout_boot,]
    Y_boot_b <- R_sqrt %*% Y_tilde_boot #RE-INDUCING DEPENDENCE
    #RD_boot[,b] <- apply(Y_boot_b,1, function(x) t(x-xbar_star) %*% invcov %*% (x-xbar_star))
    #following two lines are alternative way to get RDs
    temp <- (Y_boot_b-xbar_star_mat) %*% invcov_sqrt
    RD_boot[,b] <- rowSums(temp * temp)
  } # end loop over bootstrap samples
  RD_scaled <- scale[m] * RD_boot # scale is from Hardin & Rocke's approach (2005)
  boot_cutoffs <- t(apply(RD_scaled[notbest,],1,quantile,probs=c(1-adj_alpha)))
  
  
  # Robust Distance of the Original Data
  RD_orig <- apply(dep_data,1, function(x) t(x-xbar_star) %*% invcov %*% (x-xbar_star))
  
  
  out.RD99 <- out.RD999 <- out.RD9999 <- rep(FALSE, n)
  out.RD99[notbest] <- (RD_orig[notbest] > boot_cutoffs[1])
  out.RD999[notbest] <- (RD_orig[notbest] > boot_cutoffs[2])
  out.RD9999[notbest] <- (RD_orig[notbest] > boot_cutoffs[3])
  
  out <- data.frame(out.RD99, out.RD999, out.RD9999)
  
  names(out) <- c('0.99 quantile','0.999 quantile','0.9999 quantile')
  result <- list(outliers=out, cutoffs=boot_cutoffs)
  return(result)
}
