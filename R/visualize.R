#' Plots one or several outlyingness measures of the same type.
#'
#' @param meas A T_ x m data.frame with each column being the time-course for a 
#'  scrubbing method. Column names should identify the method as one of the following:
#'  \code{PCA_var__leverage}, \code{PCA_kurt__leverage}, \code{PCATF__leverage},
#'  \code{PCA_var__robdist}, \code{PCA_kurt__robdist},
#'  \code{PCA_var__robdist_subset}, \code{PCA_kurt__robdist_subset},
#'  \code{DVARS_DPD}, \code{DVARS_ZD}, or \code{FD}.
#' @param cuts A 1 x m data.frame with each colum being the cutoff for a 
#'  scrubbing method. Column names should be the same as those provided for \code{meas}.
#' @param name The name of the type of outlyingness measure being plotted:
#'  \code{Leverage}, \code{RobDist}, \code{RobDist Subset}, \code{FD}, \code{DVARS}.
#' @param inMCD A T_ x m data.frame indicating whether each observation was in the 
#'  MCD subset, for each method. Is only required and used if the measure is robust 
#'  distance (subset).
#' @param ... Additional arguments to ggplot: main, sub, xlab, ...
#'
#' @return A ggplot
#' 
#' @import ggplot2
#' 
#' @export
clever_plot_indiv_panel <- function(meas, cuts, name, inMCD=NULL, ...){
  args <- list(...)

  # Extra colors:
  ## "#E78AC3"
  ## "#FFD92F"
  ## "#B8B8B8" #grey, not exactly from this palette
  colors <- list(
    PCA_var__leverage = "#8DA0CB",
    PCA_kurt__leverage = "#FC8D62",
    PCATF__leverage = "#A6D854",
    PCA_var__robdist = "#8DA0CB",
    PCA_kurt__robdist = "#FC8D62",
    PCA_var__robdist_subset = "#8DA0CB",
    PCA_kurt__robdist_subset = "#FC8D62",
    DVARS_DPD = "#66C2A5",
    DVARS_ZD = "#E5C494",
    FD = "#E78AC3"
  )
  DVARS_outs_col <- "#B8B8B8"
  name_formatted <- list(
    PCA_var__leverage = "High-variance PCs",
    PCA_kurt__leverage = "High-kurtosis PCs" ,
    PCATF__leverage = "Trend-filtered PCs",
    PCA_var__robdist = "High-variance PCs",
    PCA_kurt__robdist = "High-kurtosis PCs",
    PCA_var__robdist_subset = "High-variance PCs",
    PCA_kurt__robdist_subset = "High-kurtosis PCs",
    DVARS_DPD = "DVARS Delta % D",
    DVARS_ZD = "DVARS z-score",
    FD = "Framewise Disp."
  )

  T_ <- length(meas[[1]])

  id_outs <- !is.null(cuts)

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
  any_outs <- FALSE
  if(id_outs){
    drop_line <- list()
    if(name == "DVARS"){
      DVARS_outs <- d[['DVARS_DPD']]$out & d[['DVARS_ZD']]$out
      any_outs <- any(DVARS_outs)
      if(any_outs){
        drop_line[['DVARS']] <- data.frame(
          xmin = which(DVARS_outs) - 0.5,
          xmax = which(DVARS_outs) + 0.5,
          ymin = min(c(d[['DVARS_DPD']]$meas, d[['DVARS_ZD']]$meas)),
          ymax = max(c(d[['DVARS_DPD']]$meas, d[['DVARS_ZD']]$meas))
        )
      }
    } else {
      for(n in names(meas)){
        any_outs <- any(d[[n]]$out)
        if(any_outs){
          drop_line[[n]] <- d[[n]][d[[n]]$out,]
          drop_line[[n]]$xmin <- as.numeric(rownames(drop_line[[n]])) - 0.5
          drop_line[[n]]$xmax <- as.numeric(rownames(drop_line[[n]])) + 0.5
          drop_line[[n]]$ymin <- min(0, min(data.frame(meas)))
          drop_line[[n]]$ymax <- drop_line[[n]]$meas
        }
      }
    }
  }

  # Format data as a data.frame for ggplot.
  for(i in 1:length(meas)){
    d[[names(meas)[i]]]$method = names(meas)[i]
  }
  d <- do.call(rbind, d)
  d$idx <- 1:T_

  # Get labels.
  main <- ifelse("main" %in% names(args), args$main, name)
  sub <- ifelse("sub" %in% names(args), args$sub,
    ifelse(id_outs,
           ifelse(length(drop_line) < 0, "No outliers detected.", ""),
           "(No outlier thresholding performed)"
    ))
  #xlab <- ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)")
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

  # Make ggplot.
  plt <- ggplot()
  if(any_outs){
    if(name == 'DVARS'){
      plt <- plt +
        geom_rect(data=drop_line[['DVARS']],
          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.5, fill=DVARS_outs_col)
    } else {
      for(i in 1:length(drop_line)){
        n <- names(drop_line)[i]
        plt <- plt +
          geom_rect(data=drop_line[[n]],
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.5, fill=colors[n])
        #Text label if any outlier is detected.
      }
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
    plt <- plt + geom_hline(yintercept=cuts[[n]], linetype="dashed",
      color=ifelse(length(meas)==1, 'black', colors[n]))
  }

  xticks_width <- c(1, 2, 2.5, 3, 5)
  xticks_width <- c(xticks_width, xticks_width*10, xticks_width*100, xticks_width*1000, xticks_width*10000, xticks_width*100000)
  xticks_width <- max(xticks_width[xticks_width*2 < T_*.9])
  xticks <- c(seq(from=0, to=floor(T_*.9), by=xticks_width), T_)

  plt <- plt + labs(x=xlab, y=ylab, color="Method") +
    #coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)), ylim=c(0, ylim_max*1.2)) + #fix this line
    theme(
      axis.title.x=element_blank(),
      legend.position=legend.position,
      panel.spacing.y=unit(1.5, "lines")) +
    scale_x_continuous(expand=expansion(mult = c(.01, .01)), breaks=xticks) +
    scale_y_continuous(expand=expansion(mult = c(0, .05)))

  if(name == "DVARS"){
    print("To-do: implement dual axis for DVARS...")
  }
  #if(name == "MCD Distance"){
  #  plt <- plt + facet_grid(inMCD~.)
  #}
  return(plt)
}

#' Plots the outlyingness measures from a clever result. Can support multiple panels of
#'  different outlyingness measures, but by default, it will plot only the first method.
#'
#' @param x The clever object.
#' @param methods_to_plot Either: "one" to plot only the first method; "all" to plot
#'  all the methods that were computed; or, a character vector of desired methods.
#'  Default is "one".
#' @param FD (Optional) A length T_ vector of FD values, in mm. If provided, the FD 
#'  time series plot will be added.
#' @param FD_cut (Optional) The outlier cutoff for FD. Default is 0.5 mm.
#' @param plot_title (Optional) If provided, will add a title to the plot.
#' @param ... Additional arguments to the individual plots in each panel.
#'
#' @return A ggplot
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid ggdraw
#' 
#' @method plot clever 
#' @export
plot.clever <- function(x, methods_to_plot="one", FD=NULL, FD_cut=0.5, plot_title=NULL, ...){
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

  methods <- list(
    lev = paste0(c("PCA_var", "PCA_kurt", "PCATF"), "__leverage"),
    rbd = paste(c("PCA_var", "PCA_kurt"), "robdist", sep="__"),
    rds = paste(c("PCA_var", "PCA_kurt"), "robdist_subset", sep="__"),
    DVARS = c("DVARS_DPD", "DVARS_ZD"),
    FD = c("FD")
  )
  if(is.null(methods_to_plot) | methods_to_plot=="all"){
    available_methods <- names(x$outlier_measures)
  } else if(methods_to_plot=="one"){
    available_methods <- names(x$outlier_measures)[1]
  } else {
    available_methods <- methods_to_plot
  }
  if(!is.null(FD)){ available_methods <- c(available_methods, 'FD') }
  methods <- lapply(methods, intersect, available_methods)
  methods_any <- lapply(lapply(methods, length), as.logical)

  plots <- list()
  if(methods_any$lev){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods$lev],
      cuts = x$outlier_cutoffs[methods$lev],
      name = "Leverage",
      ...
    )))
  }
  if(methods_any$rbd){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods$rbd],
      cuts = x$outlier_cutoffs[methods$rbd],
      name = "RobDist",
      inMCD = x$inMCD,
      ...
    )))
  }
  if(methods_any$rds){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods$rds],
      cuts = x$outlier_cutoffs[methods$rds],
      name = "RobDist Subset",
      inMCD = x$inMCD,
      ...
    )))
  }
  if(methods_any$DVARS){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = x$outlier_measures[methods$DVARS],
      cuts = x$outlier_cutoffs[methods$DVARS],
      name = "DVARS",
      ...
    )))
  }
  if(methods_any$FD){
    plots <- append(plots, list(clever_plot_indiv_panel(
      meas = data.frame(FD=FD),
      cuts = data.frame(FD=FD_cut),
      name = "FD",
      ...
    )))
  }

  plots[[length(plots)]] <- plots[[length(plots)]] + theme(axis.title.x=element_text()) + xlab('Index (Time Point)')

  rel_heights <- rep(1, length(plots))
  rel_heights[length(plots)] <- 1.1

  plt <- plot_grid(plotlist=plots, ncol=1, vjust=0, align="v", rel_heights=rel_heights)

  if(!is.null(plot_title)){
    plt <- plot_grid(
      ggdraw() + 
        draw_label(plot_title, fontface='bold', x=0, hjust=0) +
        theme(plot.margin = margin(0, 0, 0, 7)),
      plt,
      ncol=1,
      rel_heights=c(.15, length(plots))
    )
  }

  return(plt)
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
#' @param NA_fill Replace in-mask NA values with this. Default NA (no action).
#'
#' @return A 4D array representing the volume time series. Time is on the 4th
#'  dimension.
#'
#' @export
Matrix_to_VolumeTimeSeries <- function(mat, mask, sliced_dim = NA, NA_fill=FALSE){
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

  vts[is.na(vts)] <- NA_fill

  return(vts)
}
