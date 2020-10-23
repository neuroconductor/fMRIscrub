#' \code{"clever"} subplot
#'
#' Plots the outlyingness measures of the same type from a \code{"clever"} 
#'  result.
#'
#' @param meas A \eqn{T_ x m} data.frame with each column being the time-course for an 
#'  outlyingness measure. Column names should identify the method as one of the following:
#'  \code{PCA_var__leverage}, \code{PCA_kurt__leverage}, \code{PCATF__leverage};
#'  \code{PCA_var__robdist}, \code{PCA_kurt__robdist};
#'  \code{DVARS}, \code{DVARS_DPD}, \code{DVARS_ZD}; or \code{FD}.
#' @param cuts An \eqn{m} length vector with each value being the cutoff for a 
#'  outlyingness measure. Column names should be the same as those provided for \code{meas}.
#' @param flag An \eqn{T_ x m} data.frame with each column being the flags for an
#'  outlyingness measure.
#' @param name The name of the type of outlyingness measure being plotted:
#'  \code{Leverage}, \code{RobDist}, \code{FD}, \code{DVARS}.
#' @param robdist_info A list containing information related to the robust distance measure: inMCD, outMCD_scale, and 
#'  Fparam. Not required if robust distance-based measures are not being plottted.
#' @param ... Additional arguments to ggplot: main, sub, xlab, ...
#'
#' @return A ggplot
#' 
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' 
#' @keywords internal
clever_plot_indiv_panel <- function(meas, cuts, flag, name, robdist_info=NULL, ...){
  args <- list(...)

  xmin <- xmax <- ymin <- ymax <- idx <- method <- inMCD <- NULL

  # Extra colors:
  colors_grey <- c("#494949", "#808080", "#b8b8b8", "#cdcdcd")
  colors <- list(
    # Leverage or RobDist
    PCA_var = "#8DA0CB", # blue
    PCA_kurt = "#FC8D62", # orange
    PCATF = "#A6D854", # green
    # DVARS
    traditional = "#66C2A5", # aqua 
    DPD = "#927c5b", # dark-tan
    ZD = "#e7caa0", # light-tan
    dual = "#E5C494", # tan
    # Motion
    Motion_t1 = colors_grey[2],
    Motion_t2 = colors_grey[3],
    Motion_t3 = colors_grey[4],
    Motion_r1 = colors_grey[2],
    Motion_r2 = colors_grey[3],
    Motion_r3 = colors_grey[4],
    FD = "#E78AC3", # pink
    # GSR
    GSR = "#B8B8B8" # light-grey
  )
  colors_extra <- c(
    "#FFD92F" # yellow
  )
  name_formatted <- list(
    # Leverage or RobDist
    PCA_var = "High-variance PCs",
    PCA_kurt = "High-kurtosis PCs" ,
    PCATF = "Trend-filtered PCs",
    # DVARS
    traditional = "traditional DVARS",
    DPD = "DVARS Delta % D",
    ZD = "DVARS z-score",
    dual = "DVARS dual cutoff",
    # Motion
    Motion_t1 = "Translation RP 1",
    Motion_t2 = "Translation RP 2",
    Motion_t3 = "Translation RP 3",
    Motion_r1 = "Rotation RP 1",
    Motion_r2 = "Rotation RP 2",
    Motion_r3 = "Rotation RP 3",
    FD = "Framewise Displacement"
  )

  T_ <- nrow(meas); meas_subnames <- names(meas)

  id_outs <- !is.null(cuts)

  mcd_meas <- log_meas <- name == "RobDist"

  meas <- stack(meas)
  names(meas)[names(meas)=="values"] <- "measure"
  names(meas)[names(meas)=="ind"] <- "name"
  meas$idx <- rep(1:T_)
  # MCD
  if (mcd_meas) {
    meas$inMCD <- vector("logical", T_)
    for (ii in 1:length(meas_subnames)) {
      meas[meas$name==meas_subnames[ii],"inMCD"] <- robdist_info[meas_subnames]$inMCD
    }
  }

  if (id_outs) {
    flag <- stack(flag)
    names(flag)[names(flag)=="values"] <- "isOutlier"
    names(flag)[names(flag)=="ind"] <- "name"
    flag$idx <- rep(1:T_)
  }

  # Log values if applicable.
  if(log_meas){
    meas$measure <- log(meas$measure, base=10)
    if(id_outs){
      cuts <- log(cuts, base=10)
    }
  }

  # Get the upper y-axis limit.
  ylim_max <- ifelse(
    name=="Leverage" & (!("leverage__PCATF" %in% meas_subnames)),
    1, 
    max(max(cuts), max(meas$measure))
  )
  ylim_max <- ylim_max*1.05

  # Get the lower y-axis limit.
  ylim_min <- ifelse(name %in% c("DVARS", "GSR"), min(meas$measure), 0)

  # Check if any outliers were detected.
  if(id_outs){
    any_outs <- any(flag$isOutlier)
    drop_line <- vector("list", nrow(flag)); names(drop_line) <- names(flag)
    for (ii in 1:length(drop_line)) {
      flag_ii <- which(flag[flag$name==names(drop_line)[ii], "isOutlier"])
      if (length(flag_ii) == 0) { next }
      drop_line[[ii]] <- data.frame(
        xmin = flag_ii - 0.5, 
        xmax = flag_ii + 0.5,
        ymin = ylim_min,
        ymax = ylim_max
      )
    }
  } else {
    any_outs <- FALSE
  }

  # Get labels.
  main <- ifelse("main" %in% names(args), args$main, name)
  sub <- ifelse(
    "sub" %in% names(args), 
    args$sub,
    ifelse(
      id_outs,
      ifelse(length(drop_line) < 0, "No outliers detected.", ""),
      "(No outlier thresholding performed)"
    )
  )
  #xlab <- ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)")
  ylab <- ifelse(
    "ylab" %in% names(args), 
    args$ylab,
    ifelse(log_meas, paste0(name, " (log)"), name)
  )
  legend.position <- ifelse(
    "show.legend" %in% names(args),
    ifelse(args$show.legend, "right", "none"),
    "right"
  )

  # Make ggplot.
  plt <- ggplot()

  # Draw drop-down lines for outliers.
  if(any_outs){
    for (ii in 1:length(drop_line)) {
      if (is.null(drop_line[[ii]])) { next }
      n <- names(drop_line)[i]
      plt <- plt +
        geom_rect(
          data=drop_line[[n]],
          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.5, #fill=colors[n]
        )
      #Text label if any outlier is detected.
    }
  }

  # Outlyingness cutoff line.
  for(i in 1:length(cuts)){
    n <- names(cuts)[i]
    plt <- plt + 
      geom_hline(
        yintercept=cuts[[n]], linetype="dashed" #, color=ifelse(length(cuts)==1, "black", colors[n])
      )
  }

  # Draw data points (after drop-down lines, so they are drawn on top).
  if(mcd_meas){
    plt <- plt + 
      geom_point(data=meas, aes(x=idx, y=measure, color=name, shape=inMCD)) +
      scale_shape_manual(values=c(3, 16))
  } else {
    plt <- plt + geom_point(data=meas, aes(x=idx, y=measure, color=name))
  }
  # plt <- plt + scale_color_manual(values=colors, labels=name_formatted)

  # Use an optimal spacing between the x-ticks.
  xticks_width <- c(1, 2, 2.5, 3, 5)
  xticks_width <- c(xticks_width, xticks_width*10, xticks_width*100, 
    xticks_width*1000, xticks_width*10000, xticks_width*100000)
  xticks_width <- max(xticks_width[xticks_width*2 < T_*.9])
  xticks <- c(seq(from=0, to=floor(T_*.9), by=xticks_width), T_)

  plt <- plt + labs(y=ylab, color="Method") +
    theme_cowplot() +
    #coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)), ylim=c(0, ylim_max*1.2)) + #fix this line
    theme(
      axis.title.x=element_blank(),
      legend.position=legend.position,
      panel.spacing.y=unit(1.5, "lines")) +
    scale_x_continuous(expand=expansion(mult = c(.01, .01)), breaks=xticks) +
    scale_y_continuous(expand=expansion(mult = c(0, .01)))

  # [TO-DO]: implement dual axis for DVARS.

  return(plt)
}

#' Plot \code{"clever"}
#' 
#' Plots the outlyingness measures from a \code{"clever"} result. Can support 
#'  multiple panels of different outlyingness measures, but by default, 
#'  it will plot only the first measures.
#'
#' @param x The \code{"clever"} object.
#' @param measures "all" to plot each measure (default), or a character vector 
#'  of desired measures.
#' @param title (Optional) If provided, will add a title to the plot.
#' @param ... Additional arguments to ggplot: main, sub, xlab, ...
#'
#' @return A ggplot
#' 
#' @method plot clever 
#' @export
plot.clever <- function(x, measures="all", title=NULL, ...){
  projection_methods = x$params$projection_methods
  outlyingness_methods = x$params$outlyingness_methods
  DVARS = x$params$DVARS

  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package \"cowplot\" needed to use `plot.clever`. Please install it.", call. = FALSE)
  }

  args <- list(...)

  projection_methods_formatted <- list(
    PCA_var = "High-variance PCs",
    PCA_kurt = "High-kurtosis PCs",
    PCATF = "Trend-filtered PCs"
  )

  # Define all the subplots: measures
  measures_to_plot <- list(
    Leverage = paste0("leverage__",  c("PCA_var", "PCA_kurt", "PCATF")),
    RobDist = paste("robdist__", c("PCA_var", "PCA_kurt")),
    DVARS = c("DVARS__traditional", "DVARS__DPD", "DVARS__ZD"),
    Motion = c(paste0("motion_t", 1:3), paste0("motion_r", 1:3), "FD")
  )
  CompCor_meas <- names(x$measures)[grepl("CompCor__", names(x$measures), fixed=TRUE)]
  if (length(CompCor_meas) > 0) {
    CompCor_meas <- unique(gsub("__PC.*", "", CompCor_meas, fixed=TRUE))
    for (ii in 1:length(CompCor_meas)) {
      CompCor_meas_ii <- list(
        names(x$measures)[grepl(CompCor_meas[ii], names(x$measures), fixed=TRUE)]
      )
      names(CompCor_meas_ii) <- CompCor_meas[ii]
      measures_to_plot <- c(measures_to_plot, CompCor_meas_ii)
    }
  }
  measures_to_plot <- c(measures_to_plot, list(GSR="GSR"))

  # Define all the subplots: outliers
  outcuts_to_plot <- list(
    Leverage = paste0("leverage__",  c("PCA_var", "PCA_kurt", "PCATF")),
    RobDist = paste("robdist__", c("PCA_var", "PCA_kurt")),
    DVARS = c("DVARS__traditional", "DVARS__DPD", "DVARS__ZD"),
    Motion = "FD"
  )
  outflag_to_plot <- list(
    Leverage = paste0("leverage__",  c("PCA_var", "PCA_kurt", "PCATF")),
    RobDist = paste("robdist__", c("PCA_var", "PCA_kurt")),
    DVARS = c("DVARS__traditional", "DVARS__dual"),
    Motion = "FD"
  )

  # Remove empty or unwanted subplots
  measures_to_plot <- lapply(measures_to_plot, function(y){y[y %in% names(x$measures)]})
  outcuts_to_plot <- lapply(outcuts_to_plot, function(y){y[y %in% names(x$measures)]})
  outflag_to_plot <- lapply(
    outflag_to_plot, 
    function(y){y[y %in% names(x$measures) | (y=="DVARS__dual" & ("DVARS__DPD" %in% x$measures && "DVARS__ZD" %in% x$measures))] }
  )
  if (!("all" %in% measures)) {
    measures_to_plot <- lapply(measures_to_plot, function(y){y[y %in% measures]})
    outcuts_to_plot <- lapply(outcuts_to_plot, function(y){y[y %in% measures]})
    outflag_to_plot <- lapply(
      outflag_to_plot, 
      function(y) { y[y %in% measures | (y=="DVARS__dual" & ("DVARS__DPD" %in% measures && "DVARS__ZD" %in% measures))] }
    )
  }
  measures_to_plot <- measures_to_plot[sapply(measures_to_plot, length) > 0]

  plots <- vector("list", length(measures_to_plot))
  for (ii in 1:length(measures_to_plot)) {
    subplot_name <- names(measures_to_plot)[ii]
    plots[[ii]] <- clever_plot_indiv_panel(
      meas = x$measures[measures_to_plot[[subplot_name]]],
      cuts = x$outlier_cutoffs[outcuts_to_plot[[subplot_name]]],
      flag = x$outlier_flags[outflag_to_plot[[subplot_name]]],
      name = subplot_name,
      ...
    )
  }

  # Add x-axis label to bottom plot.
  plots[[length(plots)]] <- plots[[length(plots)]] + 
    theme(axis.title.x=element_text()) + 
    xlab(ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)"))
  rel_heights <- rep(1, length(plots))
  rel_heights[length(plots)] <- 1.1

  plt <- cowplot::plot_grid(plotlist=plots, ncol=1, vjust=0, align="v", rel_heights=rel_heights)

  # Add title if provided.
  if(!is.null(title)){
    plt <- cowplot::plot_grid(
      ggdraw() + 
        draw_label(title, fontface='bold', x=0, hjust=0) +
        theme(plot.margin = margin(0, 0, 0, 7)),
      plt,
      ncol=1,
      rel_heights=c(.15, length(plots))
    )
  }

  plt
}