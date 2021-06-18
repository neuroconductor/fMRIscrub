# #' \code{"clever"} subplot
# #'
# #' Plots the outlyingness measures of the same type from a \code{"clever"} 
# #'  result.
# #'
# #' @param meas A \eqn{T_ x m} data.frame with each column being the time-course for an 
# #'  outlyingness measure. Column names should identify the method as one of the following:
# #'  \code{PCA__leverage}, \code{PCA_kurt__leverage}, \code{PCATF__leverage};
# #'  \code{PCA__robdist}, \code{PCA_kurt__robdist};
# #'  \code{ICA__leverage}, \code{ICA_kurt__leverage};
# #'  \code{ICA__robdist}, \code{ICA_kurt__robdist};
# #'  \code{DVARS}, \code{DVARS_DPD}, \code{DVARS_ZD}; or \code{FD}.
# #' @param cuts An \eqn{m} length vector with each value being the cutoff for a 
# #'  outlyingness measure. Column names should be the same as those provided for \code{meas}.
# #' @param flag An \eqn{T_ x m} data.frame with each column being the flags for an
# #'  outlyingness measure.
# #' @param name The name of the type of outlyingness measure being plotted:
# #'  \code{PCA_leverage}, \code{PCA_robdist}, \code{ICA_leverage}, 
# #'  \code{ICA_robdist}, \code{FD}, \code{DVARS}.
# #' @param robdist_info A list containing information related to the robust 
# #'  distance measure: inMCD, outMCD_scale, and Fparam. Not used if robust 
# #'  distance-based measures are not being plottted.
# #' @param ... Additional arguments to ggplot: main, sub, xlab, ...
# #'
# #' @return A ggplot
# #' 
# #' @importFrom utils stack
# #' 
# #' @keywords internal
# clever_plot_indiv_panel <- function(meas, cuts, flag, name, robdist_info=NULL, ...){
#   colnames(meas) <- gsub("^.*__", "", colnames(meas))

#   if (!requireNamespace("cowplot", quietly = TRUE)) {
#     stop("Package \"cowplot\" needed to use `clever_plot_indiv_panel`. Please install it.", call. = FALSE)
#   }
#   if (!requireNamespace("ggplot2", quietly = TRUE)) {
#     stop("Package \"ggplot2\" needed to use `clever_plot_indiv_panel`. Please install it.", call. = FALSE)
#   }

#   args <- list(...)
#   xmin <- xmax <- ymin <- ymax <- idx <- measure <- inMCD <- NULL

#   # Extra colors:
#   colors_grey <- c("#494949", "#808080", "#b8b8b8", "#cdcdcd")
#   colors <- list(
#     # leverage or robdist
#     PCA = "#8DA0CB", # blue
#     PCA_kurt = "#FC8D62", # orange
#     PCATF = "#A6D854", # green
#     ICA = "#8DA0CB", # blue
#     ICA_kurt = "#FC8D62", # orange
#     PCA2 = "#8DA0CB", # blue
#     PCA2_kurt = "#FC8D62", # orange
#     PCATF = "#A6D854", # green
#     ICA2 = "#8DA0CB", # blue
#     ICA2_kurt = "#FC8D62", # orange
#     # DVARS
#     traditional = "#66C2A5", # aqua 
#     DPD = "#927c5b", # dark-tan
#     ZD = "#e7caa0", # light-tan
#     dual = "#E5C494", # tan
#     # motion
#     Motion_t1 = colors_grey[2],
#     Motion_t2 = colors_grey[3],
#     Motion_t3 = colors_grey[4],
#     Motion_r1 = colors_grey[2],
#     Motion_r2 = colors_grey[3],
#     Motion_r3 = colors_grey[4],
#     FD = "#E78AC3", # pink
#     # GSR
#     GSR = "#B8B8B8" # light-grey
#   )
#   colors_extra <- c(
#     "#FFD92F", # yellow
#     "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D" #Dark2
#   )
#   name_formatted <- list(
#     # leverage or robdist
#     PCA = "AAV PCs",
#     PCA_kurt = "AAV & high-kurtosis PCs" ,
#     PCATF = "AAV trend-filtered PCs",
#     ICA = "AAV ICs",
#     ICA_kurt = "AAV & high-kurtosis ICs" ,
#     PCA2 = "PESEL PCs",
#     PCA2_kurt = "PESEL & high-kurtosis PCs" ,
#     ICA2 = "PESEL ICs",
#     ICA2_kurt = "PESEL & high-kurtosis ICs" ,
#     # DVARS
#     traditional = "Traditional DVARS",
#     DPD = "DVARS Delta % D",
#     ZD = "DVARS z-score",
#     dual = "DVARS dual cutoff",
#     # motion
#     Motion_t1 = "Translation RP 1",
#     Motion_t2 = "Translation RP 2",
#     Motion_t3 = "Translation RP 3",
#     Motion_r1 = "Rotation RP 1",
#     Motion_r2 = "Rotation RP 2",
#     Motion_r3 = "Rotation RP 3",
#     FD = "Framewise Displacement"
#   )
#   ylab_formatted <- list(
#     PCA_leverage="PCA Leverage",
#     PCA_robdist="PCA Rob. Dist. (AAV)",
#     ICA_leverage="ICA Leverage",
#     ICA_robdist="ICA Rob. Dist. (AAV)",
#     PCA2_leverage="PCA Lev (PESEL)",
#     PCA2_robdist="PCA Rob. Dist. (PESEL)",
#     ICA2_leverage="ICA Lev (PESEL)",
#     ICA2_robdist="ICA Rob. Dist. (PESEL)",
#     DVARS="DVARS",
#     DVARS2="Dual DVARS",
#     motion="Motion (FD)",
#     CompCor_wm_cort="CompCor: Cortical WM",
#     CompCor_wm_cblm="CompCor: Cerebellar WM",
#     CompCor_csf="CompCor: CSF",
#     GSR="Global Signal"
#   )

#   if (!all(colnames(meas) %in% names(colors))) {
#     if (grepl("CompCor", name)) {
#       max_nPC <- max(as.numeric(gsub("PC", "", colnames(meas))))
#       if (max_nPC <= 1) {
#         colors <- c(colors, list(PC1="#000000"))
#       } else {
#         new_colors <- as.character(as.hexmode(round(seq(0, 200, length.out = max_nPC))))
#         new_colors <- paste0("#", vapply(new_colors, function(x){paste0(rep(x, 3), collapse="")}, ""))
#         names(new_colors) <- paste0("PC", seq_len(max_nPC))
#         new_colors <- as.list(new_colors)
#         colors <- c(colors, new_colors)
#       }
#     } else {
#       new_colors <- rep(colors_extra, ceiling(sum(!(colnames(meas) %in% names(colors)))/length(colors_extra)))
#       new_colors <- new_colors[seq_len(length(colors_extra))]
#       names(new_colors) <- colnames(meas)[!(colnames(meas) %in% names(colors))]
#       new_colors <- as.list(new_colors)
#       colors <- c(colors, new_colors)
#     }
#   }
#   colors <- do.call(c, colors)

#   T_ <- nrow(meas); meas_subnames <- names(meas)

#   id_outs <- !is.null(flag) && length(flag) > 0

#   mcd_meas <- log_meas <- grepl("robdist", name)

#   if (ncol(meas)==1) {
#     meas <- data.frame(
#       measure=meas[,1],
#       name=colnames(meas)
#     )
#   } else {
#     meas <- stack(meas)
#     names(meas)[names(meas)=="values"] <- "measure"
#     names(meas)[names(meas)=="ind"] <- "name"
#   }
#   meas$idx <- seq(T_)

#   # MCD
#   if (mcd_meas) {
#     meas$inMCD <- vector("logical", T_)
#     for (ii in 1:length(meas_subnames)) {
#       ii_row <- meas$name==meas_subnames[ii]
#       meas[ii_row,"inMCD"] <- robdist_info[[meas_subnames[ii]]]$inMCD
#       meas[ii_row,"measure"] <- meas[ii_row,"measure"] * 
#         ifelse(
#           meas[ii_row,"inMCD"], 
#           1, 
#           robdist_info[[meas_subnames[ii]]]$outMCD_scale
#         )
#     }
#   }

#   if (!is.null(cuts)) { names(cuts) <- gsub("^.*__", "", names(cuts)) }

#   if (id_outs) {
#     colnames(flag) <- gsub("^.*__", "", colnames(flag))
#     flag <- stack(flag)
#     names(flag)[names(flag)=="values"] <- "isOutlier"
#     names(flag)[names(flag)=="ind"] <- "name"
#     flag$idx <- rep(1:T_)
#   }

#   # Log values if applicable.
#   if(log_meas){
#     meas$measure <- log(meas$measure, base=10)
#     if(id_outs){
#       if (!is.null(cuts)) { cuts <- log(cuts, base=10) }
#     }
#   }

#   # Get the upper y-axis limit.
#   ylim_max <- ifelse(
#     grepl("Leverage", name) & (!("leverage__PCATF" %in% meas_subnames)),
#     1, 
#     ifelse(length(cuts) < 1, max(meas$measure), max(max(cuts), max(meas$measure)))
#   )
#   ylim_max <- ylim_max*1.05

#   # Get the lower y-axis limit.
#   ylim_min <- ifelse(name %in% c("DVARS", "GSR"), min(meas$measure), 
#     ifelse(log_meas, min(meas$measure), 0)
#   )

#   # Check if any outliers were detected.
#   if(id_outs){
#     any_outs <- any(flag$isOutlier)
#     drop_line <- vector("list", length(levels(flag$name))); names(drop_line) <- levels(flag$name)
#     for (ii in 1:length(drop_line)) {
#       flag_ii <- flag[flag$isOutlier & flag$name==levels(flag$name)[ii], "idx"]
#       if (length(flag_ii) < 1) { next }
#       drop_line[[ii]] <- data.frame(
#         xmin = flag_ii - 0.5, 
#         xmax = flag_ii + 0.5,
#         ymin = ylim_min,
#         ymax = ylim_max
#       )
#     }
#   } else {
#     drop_line <- NULL
#     any_outs <- FALSE
#   }

#   # Get labels.
#   main <- ifelse("main" %in% names(args), args$main, name)
#   sub <- ifelse(
#     "sub" %in% names(args), 
#     args$sub,
#     ifelse(
#       id_outs,
#       ifelse(length(drop_line) < 0, "No outliers detected.", ""),
#       "(No outlier thresholding performed)"
#     )
#   )
#   #xlab <- ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)")
#   ylab <- ifelse(name %in% names(ylab_formatted), ylab_formatted[[name]], name)
#   ylab <- ifelse(
#     "ylab" %in% names(args), 
#     args$ylab,
#     ifelse(log_meas, paste0(ylab, " (log)"), ylab)
#   )
#   legend.position <- ifelse(
#     "show.legend" %in% names(args),
#     ifelse(args$show.legend, "right", "none"),
#     "right"
#   )


#   # Make ggplot.
#   plt <- ggplot2::ggplot()

#   # Draw drop-down lines for outliers.
#   if(any_outs){
#     for (ii in seq_len(length(drop_line))) {
#       if (is.null(drop_line[[ii]])) { next }
#       n <- names(drop_line)[ii]
#       plt <- plt +
#         ggplot2::geom_rect(
#           data=drop_line[[n]],
#           ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.5, fill=colors[n]
#         )
#       #Text label if any outlier is detected.
#     }
#   }

#   # Outlyingness cutoff line.
#   for (ii in seq_len(length(cuts))) {
#     n <- names(cuts)[ii]
#     plt <- plt + 
#       ggplot2::geom_hline(
#         yintercept=cuts[n], linetype="dashed", color=ifelse(length(cuts)==1, "black", colors[n])
#       )
#   }

#   # Draw data points (after drop-down lines, so they are drawn on top).
#   if(mcd_meas){
#     plt <- plt + 
#       ggplot2::geom_point(data=meas, ggplot2::aes(x=idx, y=measure, color=name, shape=inMCD)) +
#       ggplot2::scale_shape_manual(values=c(16, 3))
  
#   # [TO DO]: Only show first 10 or so CompCor PCs, and say so in the subtitle
#   } else if (grepl("CompCor", name)) {
#     max_nPC <- max(as.numeric(gsub("PC", "", unique(meas$name))))
#     for (ii in seq(max_nPC, 1)) {
#       plt <- plt + 
#         ggplot2::geom_line(
#           data=subset(meas, name==paste0("PC", ii)), 
#           ggplot2::aes(x=idx, y=measure, group=name, color=name), size=1
#         )
#     }
#   } else if (name=="GSR") {
#     plt <- plt + 
#       ggplot2::geom_line(data=meas, ggplot2::aes(x=idx, y=measure, group=name, color=name), size=1)
#   } else {
#     plt <- plt + 
#       ggplot2::geom_point(data=meas, ggplot2::aes(x=idx, y=measure, color=name))
#   }
#   plt <- plt + ggplot2::scale_color_manual(values=colors, labels=name_formatted)

#   # Use an optimal spacing between the x-ticks.
#   xticks_width <- c(1, 2, 2.5, 3, 5)
#   xticks_width <- c(xticks_width, xticks_width*10, xticks_width*100, 
#     xticks_width*1000, xticks_width*10000, xticks_width*100000)
#   xticks_width <- max(xticks_width[xticks_width*2 < T_*.9])
#   xticks <- c(seq(from=0, to=floor(T_*.9), by=xticks_width), T_)

#   plt <- plt + ggplot2::labs(y=ylab, color="Method") +
#     cowplot::theme_cowplot() +
#     #coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)), ylim=c(0, ylim_max*1.2)) + #fix this line
#     ggplot2::theme(
#       axis.title.x=ggplot2::element_blank(),
#       legend.position=legend.position,
#       panel.spacing.y=ggplot2::unit(1.5, "lines")) +
#     ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult = c(.01, .01)), breaks=xticks) +
#     ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult = c(0, .01)))

#   # [TO-DO]: implement dual axis for DVARS.

#   return(plt)
# }

# #' Plot \code{"clever_multi"}
# #' 
# #' Plots the outlyingness measures from a \code{"clever_multi"} result. 
# #'
# #' @param x The \code{"clever_multi"} object.
# #' @param measures "all" to plot each measure (default), or a character vector 
# #'  of desired measures.
# #' @param title (Optional) If provided, will add a title to the plot.
# #' @param ... Additional arguments to ggplot: main, sub, xlab, ...
# #'
# #' @return A ggplot
# #' 
# #' @method plot clever_multi 
# #' @export
# plot.clever_multi <- function(x, measures="all", title=NULL, ...){

#   if (!requireNamespace("cowplot", quietly = TRUE)) {
#     stop("Package \"cowplot\" needed to plot the clever results. Please install it.", call. = FALSE)
#   }
#   if (!requireNamespace("ggplot2", quietly = TRUE)) {
#     stop("Package \"ggplot2\" needed to plot the clever results. Please install it.", call. = FALSE)
#   }

#   # Define all the subplots: measures
#   measures_to_plot <- list(
#     PCA_leverage = paste0("leverage__",  c("PCA", "PCA_kurt", "PCATF")),
#     PCA_robdist = paste0("robdist__", c("PCA", "PCA_kurt")),
#     PCA2_leverage = paste0("leverage__",  c("PCA2", "PCA2_kurt")),
#     PCA2_robdist = paste0("robdist__", c("PCA2", "PCA2_kurt")),
#     ICA_leverage = paste0("leverage__",  c("ICA", "ICA_kurt")),
#     ICA_robdist = paste0("robdist__", c("ICA", "ICA_kurt")),
#     ICA2_leverage = paste0("leverage__",  c("ICA2", "ICA2_kurt")),
#     ICA2_robdist = paste0("robdist__", c("ICA2", "ICA2_kurt")),
#     DVARS = "DVARS__traditional",
#     DVARS2 = c("DVARS__DPD", "DVARS__ZD"),
#     motion = c(paste0("motion_t", 1:3), paste0("motion_r", 1:3), "FD")
#   )
#   # CompCor_meas <- names(x$measures)[grepl("CompCor_", names(x$measures), fixed=TRUE)]
#   # if (length(CompCor_meas) > 0) {
#   #   CompCor_meas <- unique(gsub("__PC.*", "", CompCor_meas))
#   #   for (ii in 1:length(CompCor_meas)) {
#   #     CompCor_meas_ii <- list(
#   #       names(x$measures)[grepl(CompCor_meas[ii], names(x$measures), fixed=TRUE)]
#   #     )
#   #     names(CompCor_meas_ii) <- CompCor_meas[ii]
#   #     measures_to_plot <- c(measures_to_plot, CompCor_meas_ii)
#   #   }
#   # }
#   measures_to_plot <- c(measures_to_plot, list(GSR="GSR"))

#   # Define all the subplots: outliers
#   outcuts_to_plot <- list(
#     PCA_leverage = paste0("leverage__",  c("PCA", "PCA_kurt", "PCATF")),
#     PCA_robdist = paste0("robdist__", c("PCA", "PCA_kurt")),
#     PCA2_leverage = paste0("leverage__",  c("PCA2", "PCA2_kurt")),
#     PCA2_robdist = paste0("robdist__", c("PCA2", "PCA2_kurt")),
#     ICA_leverage = paste0("leverage__",  c("ICA", "ICA_kurt")),
#     ICA_robdist = paste0("robdist__", c("ICA", "ICA_kurt")),
#     ICA2_leverage = paste0("leverage__",  c("ICA2", "ICA2_kurt")),
#     ICA2_robdist = paste0("robdist__", c("ICA2", "ICA2_kurt")),
#     DVARS = "DVARS__traditional",
#     DVARS2 = c("DVARS__DPD", "DVARS__ZD"),
#     motion = "FD"
#   )
#   outflag_to_plot <- list(
#     PCA_leverage = paste0("leverage__",  c("PCA", "PCA_kurt", "PCATF")),
#     PCA_robdist = paste0("robdist__", c("PCA", "PCA_kurt")),
#     PCA2_leverage = paste0("leverage__",  c("PCA2", "PCA2_kurt")),
#     PCA2_robdist = paste0("robdist__", c("PCA2", "PCA2_kurt")),
#     ICA_leverage = paste0("leverage__",  c("ICA", "ICA_kurt")),
#     ICA_robdist = paste0("robdist__", c("ICA", "ICA_kurt")),
#     ICA2_leverage = paste0("leverage__",  c("ICA2", "ICA2_kurt")),
#     ICA2_robdist = paste0("robdist__", c("ICA2", "ICA2_kurt")),
#     DVARS = "DVARS__traditional",
#     DVARS2 = "DVARS__dual",
#     motion = "FD"
#   )

#   # Remove empty or unwanted subplots
#   measures_to_plot <- lapply(measures_to_plot, function(y){y[y %in% names(x$measures)]})
#   outcuts_to_plot <- lapply(outcuts_to_plot, function(y){y[y %in% names(x$measures)]})
#   outflag_to_plot <- lapply(
#     outflag_to_plot, 
#     function(y){ y[y %in% names(x$measures) | (y=="DVARS__traditional" & "DVARS" %in% names(x$measures)) | (y=="DVARS__dual" & "DVARS__DPD" %in% names(x$measures))] }
#   )
#   if (!("all" %in% measures)) {
#     measures_to_plot <- lapply(measures_to_plot, function(y){y[y %in% measures]})
#     outcuts_to_plot <- lapply(outcuts_to_plot, function(y){y[y %in% measures]})
#     outflag_to_plot <- lapply(
#       outflag_to_plot, 
#       function(y){ y[y %in% measures | (y=="DVARS__traditional" & "DVARS" %in% measures) | (y=="DVARS__dual" & "DVARS2" %in% measures)] }
#     )
#   }
#   measures_to_plot <- measures_to_plot[vapply(measures_to_plot, length, 0) > 0]

#   plots <- vector("list", length(measures_to_plot))
#   for (ii in seq_len(length(measures_to_plot))) {
#     subplot_name <- names(measures_to_plot)[ii]

#     if (subplot_name %in% names(outcuts_to_plot)) {
#       cuts_ii <- x$outlier_cutoffs[outcuts_to_plot[[subplot_name]]]
#     } else {
#       cuts_ii <- NULL
#     }
#     if (subplot_name %in% names(outflag_to_plot)) {
#       flag_ii <- x$outlier_flags[outflag_to_plot[[subplot_name]]]
#     } else {
#       flag_ii <- NULL
#     }

#     plots[[ii]] <- clever_plot_indiv_panel(
#       meas = x$measures[measures_to_plot[[subplot_name]]],
#       cuts = cuts_ii, flag = flag_ii, name = subplot_name,
#       robdist_info = x$robdist_info,
#       ...
#     )
#   }

#   # Add x-axis label to bottom plot.
#   plots[[length(plots)]] <- plots[[length(plots)]] + 
#     ggplot2::theme(axis.title.x=ggplot2::element_text()) + 
#     ggplot2::xlab(ifelse("xlab" %in% names(list(...)), list(...)$xlab, "Index (Time Point)"))
#   rel_heights <- rep(1, length(plots))
#   rel_heights[length(plots)] <- 1.1

#   if (length(plots) == 1) {
#     plt <- plots[[1]]
#     if (!is.null(title)) {
#       plt <- plt + ggplot2::ggtitle(title)
#     }
#   } else {
#     plt <- cowplot::plot_grid(plotlist=plots, ncol=1, vjust=0, align="v", rel_heights=rel_heights)

#     # Add title if provided.
#     if(!is.null(title)){
#       plt <- cowplot::plot_grid(
#         cowplot::ggdraw() + 
#           cowplot::draw_label(title, fontface='bold', x=0, hjust=0) +
#           ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7)),
#         plt,
#         ncol=1,
#         rel_heights=c(.15, length(plots))
#       )
#     }
#   }
#   plt
# }

# #' Plot \code{"clever"}
# #'
# #' @param x The \code{"clever"} object.
# #' @param title (Optional) If provided, will add a title to the plot.
# #' @param ... Additional arguments to ggplot: main, sub, xlab, ...
# #'
# #' @return A ggplot
# #' 
# #' @method plot clever 
# #' @export
# plot.clever <- function(x, title=NULL, ...){
#   plot(clever_to_multi(x)) + ggplot2::theme(legend.position = "none")
# }