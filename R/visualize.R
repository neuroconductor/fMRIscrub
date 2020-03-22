#' Visualizes the outlier distribution of a clever object.
#'
#' Prints a dot plot of the observations (x-axis) against their outlyingness (y-axis), i.e.
#'  their leverage or robust distance.
#'
#' Cutoffs for each outlier level are marked by horizontal dashed lines. Outliers are
#'  highlighted by vertical lines which extend down to the x-axis; they are colored yellow,
#'  orange and red in order of increasing outlyingness.
#'
#' If the outlyingness measure is robust distance, observations within the MCD are plotted
#'  separately from those outside the MCD. Also, the y-axes will be log10-scaled.
#'
#' @param x A clever object.
#' @param ... additional arguments to pass to \code{\link{plot}}
#' @return The clever ggplot.
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @method plot clever
#' @export
plot.clever <- function(x, ...){
	xmax = ymax = ymin = xmin = NULL
	rm(list= c("xmax", "ymax", "ymin", "xmin"))
	PCA_trend_filtering <- x$params$PCA_trend_filtering
	PCATF.formatted <- ifelse(PCA_trend_filtering, 'PCATF', 'PCA')
	choose_PCs <- x$params$choose_PCs
	choose_PCs.formatted <- switch(choose_PCs,
		kurtosis="Kurtosis",
		variance="Variance")
	method <- x$params$method
	method.formatted <- switch(method,
		leverage="Leverage",
		robdist="Robust Distance",
		robdist_subset="Robust Distance Subset")
	measure <- switch(method,
		leverage=x$leverage,
		robdist=x$robdist,
		robdist_subset=x$robdist)
	outliers <- x$outliers
	cutoffs <- x$cutoffs
	PCA_trend_filtering <- x$params$PCA_trend_filtering
	args <- list(...)

	if(is.null(outliers)){
		stop("clever did not label outliers. Run clever again with
			`id_out`=TRUE to visualize the results. ")
	}

	#Log the y-axis if the measurement is robust distance.
	log_measure <- switch(method,
		leverage=FALSE,
		robdist=TRUE,
		robdist_subset=TRUE)

	# Identify outliers and their levels of outlyingness.
	index <- 1:length(measure)
	outlier_level.num <- apply(outliers, 1, sum)  # get outlier levels as a single factor
	outlier_level.names <- c("not an outlier", colnames(outliers))
	outlier_level <- factor(outlier_level.names[outlier_level.num + 1], levels=outlier_level.names)
	d <- data.frame(index, measure, outlier_level)
	if(method %in% c("robdist","robdist_subset")){
		d$inMCD <- ifelse(x$inMCD, "In MCD", "Not In MCD")
	}

	# The plot will have lines extending downward from outliers
	#  to the x-axis.
	is_outlier <- d$outlier_level != "not an outlier"
	any_outliers <- any(is_outlier, na.rm=TRUE) #temporary na.rm
	if(any_outliers){
		# Obtain the coordinates of the outliers' lines' vertices.
		drop_line <- d[is_outlier,]
		drop_line$outlier_level <- factor(colnames(outliers)[outlier_level.num[is_outlier]], levels=colnames(outliers)) # remove 'not an outlier' level
		drop_line$xmin <- drop_line$index - .5
		drop_line$xmax <- drop_line$index + .5
		drop_line$ymin <- 0
		drop_line$ymax <- drop_line$measure
		# ggplot will draw rows from top to bottom.
		# Ordering by increasing outlier level ensures lines for the
		#  most outlying observations are drawn last, i.e. on the top layer.
		drop_line <- drop_line[order(drop_line$outlier_level),]
	}

	if(log_measure){
		# Add 1 before log transforming to ensure a positive range.
		method.formatted <- paste0("log10(", method.formatted, " + 1)")
		d$measure <- log(d$measure + 1, base = 10)
		if(any_outliers){
			drop_line$ymax <- log(drop_line$ymax + 1, base = 10)
		}
		cutoffs <- log(cutoffs + 1, base = 10)
	}

	# The lowest, middle, and highest outlier levels are colored
	#  yellow, orange, and red, respectively.
	cols <- grDevices::hsv(h=c(.1,.05,1), s=c(.6,.8,1))
	#if(any_outliers){
	#	cols <- cols[sort(unique(
	#		outlier_level.num[outlier_level.num!=0]))]
	#}

	main <- ifelse("main" %in% names(args), args$main,
		paste0("Outlier Distribution",
			ifelse(any_outliers, "", " (None Identified)")))
	sub <- ifelse("sub" %in% names(args), args$sub,
		paste0(PCATF.formatted, ", ", choose_PCs.formatted, ", ", method.formatted))
	xlab <- ifelse("xlab" %in% names(args), args$xlab, "Index (Time Point)")
	ylab <- ifelse("ylab" %in% names(args), args$ylab, method.formatted)
	legend.position <- ifelse("show.legend" %in% names(args),
		ifelse(args$show.legend, "bottom", "none"),
		"none")
	if((method=="leverage") & (!PCA_trend_filtering)){ ylim_max <- 1 }
	else { ylim_max <- max(d$measure) }


	plt <- ggplot(d, aes(x=index,y=measure, color=outlier_level))
	if(any_outliers){
		if(method=="leverage"){ nudge_y <- ylim_max * .08 }
		else { nudge_y <- ylim_max * .12 }
		plt <- plt +
			geom_rect(data=drop_line, inherit.aes=FALSE,
				aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
				fill=outlier_level), alpha=.9) +
			#geom_text(aes(label=ifelse(is_outlier, as.character(index) ,"")), size=4,
			#	nudge_y=nudge_y, check_overlap=TRUE, show.legend=FALSE)
			geom_text_repel(aes(label=ifelse(is_outlier, as.character(index) ,"")), size=4,
				nudge_y=nudge_y, show.legend=FALSE)
	}
	plt <- plt + geom_point(show.legend=FALSE) +
	scale_color_manual(values=c("grey","black","black","black")) +
	scale_fill_manual(values=cols, labels=colnames(outliers), drop=FALSE) +
	geom_hline(yintercept=cutoffs, linetype="dashed", color="gray") +
	labs(x=xlab, y=ylab, fill="Outlier Level") +
	coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)),
		ylim=c(0, ylim_max*1.2)) +
	theme_classic() +
	theme(legend.position=legend.position, panel.spacing.y=unit(1.5, "lines")) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	ggtitle(main, subtitle=sub) #+ geom_rug(sides="l", col=rgb(.5,0,0,alpha=.2))

	if(method %in% c("robdist", "robdist_subset")){
		plt <- plt + facet_grid(inMCD~.)
	}

	if("type" %in% names(args)){
		if(args$type == "n"){ return(plt) }
	}
	return(plt)
}

#'  Calculate the leverage images for each outlier that meets the
#'  \code{outlier_level} threshold, with 3 (default) being the highest/strictest
#'	and 1 being the lowest.
#'
#' @param x A clever object.
#' @param outlier_level The outlier threshold for the images, with 3 (default)
#'  being the highest/strictest. If no outliers at or above this threshold
#'  exist, no images will be made.
#'
#' @return A list of three: the mean leverage images for each outlier meeting
#'  the thresold, the top leverage images, and the indices of the top leverage
#'  images.
#'
#' @export
leverage_images <- function(x, outlier_level=3){
	svd <- x$PCs$svd
	if(is.null(svd$v)){
		stop("clever did not solve for the PC directions. Run clever
			again with `solve_directions=TRUE` to visualize the leverage images.")
	}
	N_ <- nrow(svd$v)

	outliers <- x$outliers
	if(is.null(outliers)){
		stop("clever did not label outliers. Run clever again with
			`id_out=TRUE` to visualize the results. ")
	}
	if((outlier_level < 1)|(outlier_level > 3)){
		stop("The outlier level should be 1, 2, or 3.")
	}

	lev_img_idxs <- which(outliers[,outlier_level])
	n_imgs <- length(lev_img_idxs)
	if(n_imgs == 0){
		print(paste0("clever did not find any outliers at level ",
			outlier_level, " (", colnames(outliers)[outlier_level], ")."))
		return(NULL)
	}

	lev_imgs <- list()
	lev_imgs$mean <- matrix(NA, nrow=n_imgs, ncol=N_)
	lev_imgs$top <- matrix(NA, nrow=n_imgs, ncol=N_)
	lev_imgs$top_dir <- vector(mode="numeric", length=n_imgs)
	for(i in 1:n_imgs){
		idx <- lev_img_idxs[i]
		mean_img <- svd$u[idx,] %*% t(svd$v)

		u_row <- svd$u[idx,]
		lev_imgs$mean[i,] <- u_row %*% t(svd$v)
		lev_imgs$top_dir[i] <- which.max(u_row)[1]
		lev_imgs$top[i,] <- svd$v[,lev_imgs$top_dir[i]] #Tie: use PC w/ more var.
	}

	row.names(lev_imgs$mean) <- lev_img_idxs
	row.names(lev_imgs$top) <- lev_img_idxs
	names(lev_imgs$top_dir) <- lev_img_idxs
	return(lev_imgs)
}

#'  Applies a 2D/3D mask to a matrix to get an volume time series.
#' @param mat A matrix whose rows are observations at different times, and
#'  columns are pixels/voxels.
#' @param mask A corresponding binary mask, with 1's representing regions
#'  within the area of interest and 0's representing regions to mask out.
#' @param sliced.dim If the mask is 2D, which dimension does it represent?
#'  Will default to the 3rd dimension (axial).
#'
#' @return A 4D array representing the volume time series. Time is on the 4th
#'  dimension.
#'
#' @export
Matrix_to_VolumeTimeSeries <- function(mat, mask, sliced.dim = NA){
	in.mask <- mask > 0
	t <- nrow(mat)

	if(length(dim(mask)) == 3){
		dims <- c(dim(mask), t)
	} else if(length(dim(mask)) == 2) {
		if(is.na(sliced.dim)){ sliced.dim=3 } #default to 3rd dim (axial)
		dims <- switch(sliced.dim,
									 c(1, dim(mask), t),
									 c(dim(mask)[1], 1, dim(mask)[2], t),
									 c(dim(mask), 1, t)
		)
	} else {
		stop("Not Implemented: mask must be 2D or 3D.")
	}

	vts <- array(0, dim=dims)
	for(i in 1:t){
		vts[,,,i][in.mask] <- mat[i,]
	}

	return(vts)
}
