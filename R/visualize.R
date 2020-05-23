clever_plot_indiv_panel <- function(meas, cuts, name, inMCD=NULL, ...){
  args <- list(...)
  
  # Extra colors:
  ## "#E78AC3"
  ## "#FFD92F"
  colors <- list(
    PCA_var__leverage = "#8DA0CB",
    PCA_kurt__leverage = "#FC8D62", 
    PCATF__leverage = "#A6D854", #or grey
    PCA_var__robdist = "#8DA0CB", 
    PCA_kurt__robdist = "#FC8D62",
    PCA_var__robdist_subset = "#8DA0CB",
    PCA_kurt__robdist_subset = "#FC8D62", 
    DVARS_DPD = "#66C2A5",
    DVARS_ZD = "#E5C494"
  )
  name_formatted <- list(
    PCA_var__leverage = "High-variance PCs",
    PCA_kurt__leverage = "High-kurtosis PCs" ,
    PCATF__leverage = "Trend-filtered PCs", #or grey
    PCA_var__robdist = "High-variance PCs",
    PCA_kurt__robdist = "High-kurtosis PCs",
    PCA_var__robdist_subset = "High-variance PCs",
    PCA_kurt__robdist_subset = "High-kurtosis PCs",
    DVARS_DPD = "DVARS Delta % D",
    DVARS_ZD = "DVARS z-score"
  )
  
  T_ <- length(meas[[1]])
  
  id_outs <- !is.null(cuts)
  
  #xmax = ymax = ymin = xmin = NULL
  #rm(list= c("xmax", "ymax", "ymin", "xmin"))

  if(name == "DVARS"){ 
    DVARS_names <- list(DVARS_DPD="DVARS Delta Pct. D", 
                        DVARS_ZD="DVARS z-score")
  }  
  
  mcd_meas <- c("RobDist", "RobDist Subset")
  log_meas <- name %in% mcd_meas
  if(log_meas){
    for(i in 1:length(meas)){
      n <- names(meas)[i]
      meas[[n]] <- log(meas[[n]], base = 10)
    }
    if(id_outs){ 
      for(cut in cuts){ cut <- log(cut+1, base=10) }
    }
  }
  
  # For each measure, collect relevant information into a dataframe.
  d <- list()
  for(i in 1:length(meas)){
    n <- names(meas)[i]
    d[[n]] <- data.frame(meas=meas[[n]])
    
    if(id_outs){
      d[[n]]$out <- meas[[n]] > cuts[[n]]
    }
    
    if(name %in% mcd_meas){
      d[[n]]$inMCD <- ifelse(inMCD[[n]], "In MCD", "Not In MCD")
    }
  }
  
  # For each measure, collect outlier information if any exist.
  if(id_outs){
    drop_line <- list()
     for(n in names(meas)){
      any_outs <- any(d[[n]]$out)
      if(any_outs){
        drop_line[[n]] <- d[[n]][d[[n]]$out,]
        drop_line[[n]]$xmin <- as.numeric(rownames(drop_line[[n]])) - 0.5
        drop_line[[n]]$xmax <- as.numeric(rownames(drop_line[[n]])) + 0.5
        drop_line[[n]]$ymin <- 0
        drop_line[[n]]$ymax <- drop_line[[n]]$meas
      }
    }
  }
  
  for(i in 1:length(meas)){
    d[[names(meas)[i]]]$method = names(meas)[i]
  }
  d <- do.call(rbind, d)
  d$idx <- 1:T_
  
  main <- ifelse("main" %in% names(args), args$main, name)
  sub <- ifelse("sub" %in% names(args), args$sub,
    ifelse(id_outs,
           ifelse(length(drop_line) < 0, "No outliers detected.", ""),
           "(No outlier thresholding performed)"
    ))
  xlab <- ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)")
  ylab <- ifelse("ylab" %in% names(args), args$ylab, 
    ifelse(log_meas, paste0("log(", name, " + 1)"), name))
  legend.position <- ifelse("show.legend" %in% names(args),
    ifelse(args$show.legend, "right", "none"),
    "none")
  if((name=="Leverage") & (!("PCATF__lev" %in% names(meas)))){ 
    ylim_max <- 1 
  } else { 
    ylim_max <- max(d$measure) 
  }
  
  plt <- ggplot()
  if(id_outs){
    for(i in 1:length(drop_line)){
      n <- names(drop_line)[i]
      plt <- plt +
        geom_rect(data=drop_line[[n]],
        aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.5, fill=colors[n])
      #Text label if any outlier is detected.
    }
  }
  if(name %in% mcd_meas){
    plt <- plt + geom_point(data=d, aes(x=idx, y=meas, color=method, shape=inMCD)) + 
      scale_shape_manual(values=c(3, 16))
    #plt <- plt + geom_point(data=d, aes(x=idx, y=meas, color=method))
  } else {
    plt <- plt + geom_point(data=d, aes(x=idx, y=meas, color=method))
  }
  plt <- plt +
    scale_color_manual(values=colors, labels=name_formatted)
    #scale_fill_manual(values=colors)
  for(i in 1:length(meas)){
    n <- names(meas)[i]
    plt <- plt + geom_hline(yintercept=cuts[[n]], linetype="dashed", color=colors[n])
  }
  plt <- plt + labs(x=xlab, y=ylab, color="Method") +
    #coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)), ylim=c(0, ylim_max*1.2)) + #fix this line
    theme(legend.position=legend.position, panel.spacing.y=unit(1.5, "lines")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  if(name == "DVARS"){
    print("Need to implement dual axis...")
    print("Minimum is not zero...")
  }
  #if(name == "MCD Distance"){
  #  plt <- plt + facet_grid(inMCD~.)
  #}
  return(plt)
}

plot.clever <- function(x, ...){
  projection_methods = x$params$projection_methods
  outlyingness_methods = x$params$outlyingness_methods
  DVARS = x$params$DVARS
  PCATF_kwargs = x$params$PCATF_kwargs
  kurt_quantile = x$params$kurt_quantile
  kurt_detrend = x$params$kurt_detrend
  #id_outs = x$params$id_outs
  #lev_cutoff = x$params$lev_cutoff
  #MCD_cutoff = x$params$MCD_cutoff
  lev_images = x$params$lev_images

  projection_methods_formatted <- list(
    PCA_var = "High-variance PCs",
    PCA_kurt = "High-kurtosis PCs",
    PCATF = "PCA Trend Filtering"
  )

  methods_lev <- paste0(c("PCA_var", "PCA_kurt", "PCATF"), "__leverage")
  methods_lev <- names(x$outlier_measures[names(x$outlier_measures) %in% methods_lev])
  any_lev <- length(methods_lev) > 0
  methods_rbd <- paste(c("PCA_var", "PCA_kurt"), "robdist", sep="__")
  any_rbd <- length(methods_rbd) > 0
  methods_rds <- paste(c("PCA_var", "PCA_kurt"), "robdist_subset", sep="__")
  any_rds <- length(methods_rds) > 0
  methods_DVARS <- c("DVARS_DPD", "DVARS_ZD")
  methods_DVARS <- names(x$outlier_measures[names(x$outlier_measures) %in% methods_DVARS])
  any_DVARS <- length(methods_DVARS) > 0

  plots <- list()
  if(any_lev){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods_lev],
      cuts = x$outlier_cutoffs[methods_lev],
      name = "Leverage",
      ...
    )))
  }
  if(any_rbd){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods_rbd],
      cuts = x$outlier_cutoffs[methods_rbd],
      name = "Robust Distance",
      inMCD = x$inMCD,
      ...
    )))
  }
  if(any_rds){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods_rds],
      cuts = x$outlier_cutoffs[methods_rds],
      name = "Robust Distance (Subset)",
      inMCD = x$inMCD,
      ...
    )))
  }
  if(any_DVARS){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods_DVARS],
      cuts = x$outlier_cutoffs[methods_DVARS],
      name = "DVARS",
      ...
    )))
  }
  
  plot_grid(plotlist=plots, ncol=1, align="v")
}

#' Calculate the leverage images for each outlier that meets the
#'  \code{outlier_level} threshold, with 3 (default) being the highest/strictest
#'  and 1 being the lowest.
#'
#' @param X_svd A singular value decomposition (list with u, v, and d).
#' @param timepoints The times for which to compute leverage images (rows of U).
#' @param const_mask Mask for the voxels that were removed.
#'
#' @return A list of three: the mean leverage image for each outlier meeting
#'  the thresold, the top leverage image for each outlier, and the indices of
#'  the top leverage images.
#'
#' @export
get_leverage_images <- function(X_svd, timepoints, const_mask=NULL){
  if(is.null(const_mask)){ const_mask = rep(FALSE, nrow(X_svd$v)) }
  N_ <- length(const_mask)
  n_imgs <- length(timepoints)
  lev_imgs <- list(mean=NULL, top=NULL, top_dir=NULL)
  if(n_imgs > 0){
    lev_imgs <- list()
    lev_imgs$mean <- matrix(NA, nrow=n_imgs, ncol=N_)
    lev_imgs$top <- matrix(NA, nrow=n_imgs, ncol=N_)
    lev_imgs$top_dir <- vector(mode="numeric", length=n_imgs)
    for(i in 1:n_imgs){
      idx <- timepoints[i]
      mean_img <- X_svd$u[idx,] %*% t(X_svd$v)
      u_row <- X_svd$u[idx,]
      lev_imgs$mean[i,!const_mask] <- u_row %*% t(X_svd$v)
      lev_imgs$top_dir[i] <- which.max(u_row)[1]
      lev_imgs$top[i,!const_mask] <- X_svd$v[,lev_imgs$top_dir[i]] #Tie: use PC w/ more var.
    }
    row.names(lev_imgs$mean) <- timepoints
    row.names(lev_imgs$top) <- timepoints
    names(lev_imgs$top_dir) <- timepoints
  }
  return(lev_imgs)
}

#' Applies a 2D/3D mask to a matrix to get a 3D/4D volume time series.
#' @param mat A matrix whose rows are observations at different times, and
#'  columns are pixels/voxels.
#' @param mask A corresponding binary mask, with 1's representing regions
#'  within the area of interest and 0's representing regions to mask out.
#' @param sliced_dim If the mask is 2D, which dimension does it represent?
#'  Will default to the 3rd dimension (axial).
#'
#' @return A 4D array representing the volume time series. Time is on the 4th
#'  dimension.
#'
#' @export
Matrix_to_VolumeTimeSeries <- function(mat, mask, sliced_dim = NA){
  in_mask <- mask > 0
  t <- nrow(mat)

  if(length(dim(mask)) == 3){
    dims <- c(dim(mask), t)
  } else if(length(dim(mask)) == 2) {
    if(is.na(sliced_dim)){ sliced_dim=3 } #default to 3rd dim (axial)
    dims <- switch(sliced_dim,
                   c(1, dim(mask), t),
                   c(dim(mask)[1], 1, dim(mask)[2], t),
                   c(dim(mask), 1, t)
    )
  } else {
    stop("Not Implemented: mask must be 2D or 3D.")
  }

  vts <- array(0, dim=dims)
  for(i in 1:t){
    vts[,,,i][in_mask] <- mat[i,]
  }

  return(vts)
}
