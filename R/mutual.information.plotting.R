#' Create mutual information violinplot
#'
#' Designed for visualising the mutual information matrix returned by calculate.mutual.information.pairwise.data.types.
#' @param MI matrix of mutual information as output by calculate.mutual.information.pairwise.data.types.
#' The column names should have the format <data.type.A>_<data.type.B> i.e. indicate the two variables for which mutual information was calculated delimted by '_'.
#' Rows correspond to the entities from which mutual information was calculated. i.e. genes.bounded between 0-1
#' @param log.scale plot mutual information on a log scale? Assumes mutual information is normalized (0)
#' @param ylabel label to use for y-axis . Defaults to 'Mutual Information').
#' @param yaxis.lab see ?create.violinplot
#' @param yaxis.cex see ?create.violinplot
#' @param y.min.log minimum value to display if using a log10 scale. Defaults to 0.01
#' @param style BoutrosLab.plotting style, see ?create.violinplot for more. Defaults to 'Nature'.
#' @return a violinplot displaying the distribution of observed mutual information across genes between pairs of data types.
#' @export
create.mutual.information.violinplot <- function(
	MI,
	log.scale,
	ylabel = 'Mutual Information',
	yaxis.cex = 0.55,
	y.min.log = 0.01,
	style = 'Nature',
	yaxis.lab = NULL,
	yat = NULL,
	ylimits = NULL
	) {

	MI.long <- tidyr::pivot_longer(
		data.frame(MI),
		cols = colnames(MI),
		names_to = 'data.type',
		values_to = 'MI',
		names_transform = list(data.type = ~ readr::parse_factor(.x)) #cast data.type into factor
		);
	custom.yaxis.lab <- !is.null(yaxis.lab);
	custom.yat <- !is.null(yat);

	if (log.scale) {
		if (!custom.yaxis.lab) {
			yaxis.lab <- seq(0, signif(max(MI.long$MI, na.rm = T), 1) + .1, .1);
			}
		if (!custom.yat) {
			yat <- log10(yaxis.lab);
			}
		if (!custom.yaxis.lab) {
			yaxis.lab <- c(bquote('' <= .(y.min.log)), as.expression(yaxis.lab[-1]));
			}
		if (!custom.yat) {
			yat[1] <- log10(y.min.log);
			}
		MI.long$MI[MI.long$MI <= y.min.log] <- y.min.log;
		MI[MI <= y.min.log] <- y.min.log;
		MI.medians <- apply(MI, 2, median, na.rm = TRUE);
		MI.long$MI <- log10(MI.long$MI);
		MI.medians <- log10(MI.medians);
		} else {
		MI.medians <- apply(MI, 2, median, na.rm = TRUE);
		if (!custom.yat) {
			yat <- seq(0, signif(max(MI.long$MI, na.rm = T), 1) + .1, .1);
			}
		if (!custom.yaxis.lab) {
			yaxis.lab <- yat;
			}
		}
	if (is.null(ylimits)) {
		ylimits <- range(yat);
		}

	BoutrosLab.plotting.general::create.violinplot(
		formula = MI ~ data.type,
		data = MI.long,
		style = style,
		xlab.label = NULL,
		ylab.label = bquote(bold(.(ylabel))),
		xlab.cex = 1,
		ylab.cex = 1,
		xaxis.lab = rep('', ncol(MI)),
		xaxis.tck = 0,
		yaxis.tck = 0,
		yaxis.cex = yaxis.cex,
		extra.points = list(MI.medians),
		extra.points.pch = 21,
		extra.points.col = c('white'),
		extra.points.cex = 1,
		ylimits = ylimits,
		yat = yat,
		yaxis.lab = yaxis.lab
		);
	}

#' Create a data type pair heatmap
#'
#' Create a heatmap providing a visual representation of the pariwise data type comparisons for which mutual information was calculated with calculate.mutual.information.pairwise.data.types, along with an accompanying legend.
#' @param data.type.pairs column names from the matrix of mutual information as output by calculate.mutual.information.pairwise.data.types.
#' This should be a character vector of the format <data.type.A>_<data.type.B> i.e. indicate the two variables for which mutual information was calculated delimted by '_'.
#' The vector should be ordered in a way consistent with other visualisations of the mutual information matrix e.g. by median mutual information of data type pairs.
#' @param data.type.colours a named character vector list mapping data types to desired colours. Names should match those in data.type.pairs
#' @return a list with elements data.types.map and data.types.legend.
#' data.types.map is a heatmap providing a visual representation of the pariwise data type comparisons for which mutual information was calculated with calculate.mutual.information.pairwise.data.types.
#' data.types.legend is a legend visualising the data type to colour mapping in data.type.colours.
#' @export
create.data.type.colourlegend.heatmap <- function(
	data.type.pairs,
	data.type.colours,
	data.type.legend.labels,
	yaxis.lab = c('Data type 1', 'Data type 2'),
	col.colour = 'white',
	grid.col = TRUE
	) {

	#split data type pair strings by '_' delimiter.
	data.types.colour.matrix <- sapply(
		strsplit(data.type.pairs, '_'),
		function(data.type)
			data.type.colours[rev(data.type)] #reverse order for plotting, the first element will be the bottom row of the heatmap
		);

	data.types.map <- BoutrosLab.plotting.general::create.heatmap(
		x = data.types.colour.matrix,
		clustering.method = 'none',
		yaxis.tck = 0,
		xaxis.tck = 0,
		same.as.matrix = ncol(data.types.colour.matrix) > 1,
		print.colour.key = FALSE,
		input.colours = TRUE,
		axes.lwd = 1,
		yaxis.lab = yaxis.lab,
		yaxis.cex = 0.75,
		col.colour = col.colour,
		grid.col = grid.col
		);

	data.types.legend <- list(
		colours = rev(data.type.colours),
		labels = rev(data.type.legend.labels),
		title = bquote(bold(underline('Data type'))),
		cex = .7
		);

	return(list(
		plot = data.types.map,
		legend = data.types.legend
		));
	}

#' Create heatmap displaying summary statistic of interest
#'
#' Create heatmap displaying summary statistic of interest for each pairwise data type comparison.
#' e.g. the total number of genes analysed in each pairwise data type comparison, or the number of those genes having significant mutual information can be displayed.
#' @param summary.statistic vector containing summary statistic relating to each pairwise data type compraison.
#' @param colour.scheme vector of length two mapping values of summary statistic.
#' First element corresponds to colour used for 0 value, second element corresponds to colour used for max(summary.statistic) scheme
#' @return a list with elements summary.statistic.heatmap summary.statistic.legend
#' summary.statistic.heatmap is a heatmap displaying summary.statistic.
#' summary.statistic.legend is a legend visualising the mapping of summary.statistic colour.scheme.
#' @export
create.mutual.information.summary.statistic.heatmap <- function(
	summary.statistic,
	legend.name,
	yaxis.lab,
	colour.scheme,
	col.colour = 'white',
	grid.col = TRUE,
	max.value = NULL,
	total.colours = 50
	) {
	summary.statistic.matrix <- t(data.frame(summary.statistic));
	#handle NA cases (e.g. if `num.genes` == 0)
	summary.statistics.text <- ifelse(is.na(summary.statistic.matrix), '', summary.statistic.matrix);
	if (is.null(max.value)) {
	  max.value <- max(summary.statistic, na.rm = TRUE)
	}
	min.value <- 0
	at <- seq(from = min.value, to = max.value, length.out = total.colours);
	if (ncol(summary.statistic.matrix) > 1 & length(unique(as.numeric(summary.statistic.matrix))) > 1) {
  	summary.statistic.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  		x = summary.statistic.matrix,
  		cell.text = summary.statistics.text,
  		text.col = 'black',
  		text.cex = 0.65,
  		row.pos = rep(1, length(summary.statistic)),
  		col.pos = 1:length(summary.statistic),
  		col.colour = col.colour,
  		grid.col = grid.col,
  		cluster.dimensions = 'none',
  		yaxis.tck = 0,
  		xaxis.tck = 0,
  		yaxis.lab = yaxis.lab,
  		yat = 1,
  		yaxis.cex = 0.75,
  		same.as.matrix = TRUE,
  		print.colour.key = FALSE,
  		colour.centering.value = 0,
  		colour.scheme = colour.scheme,
  		at = at,
  		total.colours = total.colours,
  		axes.lwd = 1
  		);
	  } else {
	    summary.statistic.matrix <- matrix(
	      colour.scheme[2],
	      nrow = nrow(summary.statistic.matrix),
	      ncol = ncol(summary.statistic.matrix)
	      )
	    summary.statistic.heatmap <- BoutrosLab.plotting.general::create.heatmap(
	      x = summary.statistic.matrix,
	      cell.text = summary.statistics.text,
	      input.colours = TRUE,
	      text.col = 'black',
	      text.cex = 0.65,
	      row.pos = rep(1, length(summary.statistic)),
	      col.pos = 1:length(summary.statistic),
	      col.colour = col.colour,
	      grid.col = grid.col,
	      cluster.dimensions = 'none',
	      clustering.method = 'none',
	      yaxis.tck = 0,
	      xaxis.tck = 0,
	      yaxis.lab = yaxis.lab,
	      yat = 1,
	      yaxis.cex = 0.75,
	      same.as.matrix = TRUE,
	      print.colour.key = FALSE,
	      axes.lwd = 1
	    );
    }

	#determine label breaks - make three break points
	legend.label.max <- signif(max.value, 1);
	legend.labels <- seq(min.value, max.value, legend.label.max / 2);
	if (2 == length(legend.labels)) {
		digits <- 10^(ceiling(log10(abs(legend.label.max))) - 1);
		legend.label.max <- floor(legend.label.max / digits) * digits;
		legend.labels <- seq(min.value, legend.label.max, legend.label.max / 2);
		}
	legend.at <- legend.labels / max.value * 100;
	summary.statistic.legend <- list(
		colours = colour.scheme,
		labels = legend.labels,
		title = bquote(bold(underline(.(legend.name)))),
		title.fontface = 'plain',
		continuous = TRUE,
		tck = 0,
		tck.number = 2,
		cex = .5,
		at = legend.at,
		height = 1.75
		);

	return(list(
		plot = summary.statistic.heatmap,
		legend = summary.statistic.legend
		));
	}

#' Assemble mutual information mutlipanelplot
#'
#' Combine all of the individual plots and legend into a single plot.
#' @param MI.violinplot mutual information violinplot at returned by create.mutual.information.violinplot
#' @param num.genes.heatmap a heatmap displaying numbers of genes analysed in each pairwise data type comparison. As returned by create.mutual.information.summary.statistic.heatmap.
#' @param num.sig.genes.heatmap a heatmap displaying numbers of genes with significant mutual information in each pairwise data type comparison. As returned by create.mutual.information.summary.statistic.heatmap.
#' @param plot.objects.heights relative heights of each plot i.e. in decimal scale where the sum is equal to one.
#' @param data.types.legend a legend visualising the data type to colour mapping in data.type.colours. As returned by create.data.type.colourlegend.heatmap.
#' @param num.genes.legend a legend mapping numbers of genes displayed in num.genes.heatmap to colour scheme. As returned by create.mutual.information.summary.statistic.heatmap.
#' @param num.sig.genes.legend a legend mapping numbers of genes displayed in num.sig.genes.heatmap to colour scheme. As returned by create.mutual.information.summary.statistic.heatmap.
#' @param return.plot should the plot object be returned. Defaults to FALSE.
#' @param plot.legend whether to add a legend to the plot detailing colour scheme for MI summary statistics heatmap. Defaults to TRUE.
#' @param filename path to write plot to. If NULL no image is saved. Defaults to NULL.
#' @param plot.width if filename is not NULL, the width of the image to be save in inches. See ?create.multipanelplot.
#' @param plot.height if filename is not NULL, the height of the image to be save in inches. See ?create.multipanelplot.
#' @param y.spacing vertical spacing between each plot. See ?create.multipanelplot
#' @param style 'BoutrosLab' or 'Nature'. Defaults to 'Nature'. See ?create.multipanelplot.
#' @return a multipanelplot visualising the input plots and legends.
#' @export
create.mutual.information.multipanelplot <- function(
	MI.violinplot,
	data.types.heatmap,
	num.genes.heatmap,
	num.sig.genes.heatmap,
	perc.sig.genes.heatmap,
	plot.objects.heights = c(0.56, 0.14, 0.1, 0.1, 0.1),
	data.types.legend,
	num.genes.legend,
	num.sig.genes.legend,
	perc.sig.genes.legend,
	return.plot = FALSE,
	plot.legend = TRUE,
	filename = NULL,
	plot.width = NULL,
	plot.height = NULL,
	y.spacing = -0.5,
	style = 'Nature',
	main = main,
	main.cex = main.cex
	) {

	#combine plots
	combined.plots <- list(
		MI.violinplot,
		data.types.heatmap,
		num.genes.heatmap,
		num.sig.genes.heatmap,
		perc.sig.genes.heatmap
		);
	#combine legends
	if (plot.legend) {
  	cov.legend.grob <- BoutrosLab.plotting.general::legend.grob(
  		legends = c(
  			list(
  				legend = data.types.legend,
  				legend = num.genes.legend,
  				legend = num.sig.genes.legend,
  				legend = perc.sig.genes.legend
  				)
  			),
  		label.cex = 0.5,
  		title.cex = 0.6,
  		title.just = 'left',
  		title.fontface = 'plain',
  		size = 1.5
  		);
	} else {
	  cov.legend.grob <- NULL
  }

	#final composite plot
	MI.composite.plot <- BoutrosLab.plotting.general::create.multipanelplot(
		combined.plots,
		filename = filename,
		layout.width = 1,
		layout.height = length(plot.objects.heights),
		size.units = 'in',
		resolution = 500,
		style = style,
		main = main,
		main.cex = main.cex,
		plot.objects.heights = plot.objects.heights,
		plot.objects.widths = 1,
		y.spacing = y.spacing,
		x.spacing = -0.35,
		layout.skip = rep(FALSE, length(plot.objects.heights)),
		width = plot.width,
		height = plot.height,
		bottom.padding = -2,
		top.padding = -1,
		left.legend.padding = 1,
		right.padding = 1,
		legend = list(
			right = list(
				fun = cov.legend.grob
				)
			),
		);
	if (return.plot) {
		return(MI.composite.plot);
		}
	}

#' Create a summary plot from a mutual information analysis
#'
#' Creates four plot elements and optionally assembles them into a multipanelplot or returns them individually.
#' @param MI matrix of mutual information as output by calculate.mutual.information.pairwise.data.types.
#' The column names should have the format <data.type.A>_<data.type.B> i.e. indicate the two variables for which mutual information was calculated delimted by '_'.
#' Rows correspond to the entities from which mutual information was calculated. i.e. genes.
#' @param order.by.median Should the data corresponding to each pair of data types analysed be displayed in order of median mutual information (highest to lowest). Defaults to TRUE.
#' @param log.scale plot mutual information on a log scale? Assumes mutual information is normalized (0)
#' @param data.type.colours a named character vector list mapping data types to desired colours. Names should match data types in colnames of `MI`.
#' @param summary.statistics matrix of MI summary statistics as returned by `calculate.mutual.information.pairwise.data.types`. Assumed to have row.names corresponding to colnames of `MI`.
#' Should have columns named 'num.genes' and 'num.sig.genes'.
#' @param return.plot should the plot object be returned. If return.plot.elements = TRUE, this argument is ignored. Defaults to FALSE.
#' @param return.plot.elements should individual plot elements be returned instead of assembling into a multipanelplot. If TRUE, a list with elements 'plots' and 'legends' is returned. Defaults to FALSE.
#' @param plot.legend whether to add a legend to the plot detailing colour scheme for MI summary statistics heatmap. Defaults to TRUE.
#' @param filename path to write plot to. If NULL no image is saved. If return.plot.elements = TRUE, this argument is ignored. Defaults to NULL.
#' @param plot.width if filename is not NULL, the width of the image to be saved in inches. See ?create.multipanelplot. Defaults to ncol(MI).
#' @param plot.height if filename is not NULL, the height of the image to be saved in inches. See ?create.multipanelplot. Defaults to 5.
#' @param style 'BoutrosLab' or 'Nature'. Defaults to 'Nature'. See ?create.multipanelplot.
#' @param violinplot.yaxis.cex yaxis.cex for mutual information violinplot.
#' @param violinplot.ylabel ylabel for mutual information violinplot.
#' @param violinplot.y.min.log minimum value to display in violinplot if using a log10 scale. Defaults to 0.01
#' @param violinplot.yaxis.lab yaix.lab for mutual information violinplot.
#' @param violinplot.yat yat for for mutual information violinplot.
#' @param max.value.summary.statistics a named list with elements "num.genes", "num.sig.genes", "perc.sig.genes". Used to set the maximum value of the colour scale for MI summary statistic heatmaps.
#' This parameter can be used to adjust the colour scheme and make MI plots from different datasets comparable.
#' @return Plots consist of:
#' A violinplot displaying the distribution of observed mutual information across genes between pairs of data types.
#' A heatmap providing a visual representation of the pairwise data type comparisons for which mutual information was calculated.
#' A heatmap displaying numbers of genes analysed in each pairwise data type comparison.
#' A heatmap displaying numbers of genes with significant mutual information in each pairwise data type comparison.
#' If return.plot.elements is TRUE, individual plot elements are returned.
#' If FALSE, plot elements are assembled into a multipanelplot.
#' If return.plot is TRUE, the multipanelplot is returned.
#' If return.plot.elements is FALSE and filename is not NULL, the multipanelplot is saved to a file.
#' @export
#' @examples
#' # See example code in ?calculate.mutual.information.pairwise.data.types
#'
create.mutual.information.summary.plot <- function(
	MI,
	order.by.median = TRUE,
	log.scale = TRUE,
	data.type.colours,
	data.type.legend.labels = names(data.type.colours),
	summary.statistics,
	return.plot = FALSE,
	return.plot.elements = FALSE,
	plot.legend = TRUE,
	filename = NULL,
	plot.width = ncol(MI),
	plot.height = 5,
	style = 'Nature',
	main = '',
	main.cex = 1,
	violinplot.yaxis.cex = 0.55,
	violinplot.ylabel = 'Mutual Information',
	violinplot.y.min.log = 0.01,
	violinplot.yaxis.lab = NULL,
	violinplot.yat = NULL,
	violinplot.ylimits = NULL,
	y.spacing = -0.5,
	col.colour = 'white',
	grid.col = TRUE,
	max.value.summary.statistics = list(
	  num.genes = max(summary.statistics$num.genes, na.rm = TRUE),
	  num.sig.genes = max(summary.statistics$num.sig.genes, na.rm = TRUE),
	  perc.sig.genes = 100 * max(summary.statistics$num.sig.genes / summary.statistics$num.genes, na.rm = TRUE)
	  )
	) {

	if (class(MI) != 'matrix') {
		MI <- as.matrix(MI) # functions assume MI is matrix e.g. use of colnames() rather than names()
		}

	if (order.by.median) {
		MI <- MI[, order(apply(MI, 2, median, na.rm = TRUE), decreasing = TRUE), drop = FALSE];
		summary.statistics <- summary.statistics[colnames(MI), ];
		}

	MI.violinplot <- create.mutual.information.violinplot(
		MI,
		log.scale = log.scale,
		yaxis.cex = violinplot.yaxis.cex,
		ylabel = violinplot.ylabel,
		y.min.log = violinplot.y.min.log,
		yaxis.lab = violinplot.yaxis.lab,
		yat = violinplot.yat,
		ylimits = violinplot.ylimits
		);

	data.types.plot.objects <- create.data.type.colourlegend.heatmap(
		data.type.pairs = colnames(MI),
		data.type.colours = data.type.colours,
		data.type.legend.labels = data.type.legend.labels,
		grid.col = grid.col,
		col.colour = col.colour
		);
	num.genes.plot.objects <- create.mutual.information.summary.statistic.heatmap(
		summary.statistic = summary.statistics[, 'num.genes'],
		legend.name = 'Analyzed Genes',
		yaxis.lab = 'n-Analyzed',
		colour.scheme = c('white', 'orange'),
		grid.col = grid.col,
		col.colour = col.colour,
		max.value = max.value.summary.statistics$num.genes
		);
	num.sig.genes.plot.objects <- create.mutual.information.summary.statistic.heatmap(
		summary.statistic = summary.statistics[, 'num.sig.genes'],
		legend.name = 'Significant Genes',
		yaxis.lab = 'n-Significant',
		colour.scheme = c('white', 'pink2'),
		grid.col = grid.col,
		col.colour = col.colour,
		max.value = max.value.summary.statistics$num.sig.genes
		);
	perc.sig.genes <- signif(100 * (summary.statistics[, 'num.sig.genes'] / summary.statistics[, 'num.genes']), digits = 2);
	perc.sig.genes.plot.objects <- create.mutual.information.summary.statistic.heatmap(
		summary.statistic = perc.sig.genes,
		legend.name = '% Significant Genes',
		yaxis.lab = '%Significant',
		colour.scheme = c('white', 'steelblue1'),
		grid.col = grid.col,
		col.colour = col.colour,
		max.value = max.value.summary.statistics$perc.sig.genes
		);

	if (!return.plot.elements) {

		MI.summary.plot <- create.mutual.information.multipanelplot(
			MI.violinplot = MI.violinplot,
			data.types.heatmap = data.types.plot.objects$plot,
			num.genes.heatmap = num.genes.plot.objects$plot,
			num.sig.genes.heatmap = num.sig.genes.plot.objects$plot,
			perc.sig.genes.heatmap = perc.sig.genes.plot.objects$plot,
			data.types.legend = data.types.plot.objects$legend,
			num.genes.legend = num.genes.plot.objects$legend,
			num.sig.genes.legend = num.sig.genes.plot.objects$legend,
			perc.sig.genes.legend = perc.sig.genes.plot.objects$legend,
			return.plot = return.plot,
			plot.legend = plot.legend,
			filename = filename,
			plot.width = plot.width,
			plot.height = plot.height,
			style = 'Nature',
			y.spacing = y.spacing,
			main = main,
			main.cex = main.cex
			);
		if (return.plot) {
			return(MI.summary.plot);
			}
		} else {
			return(list(
				plots = list(
					MI.violinplot = MI.violinplot,
					data.types.heatmap = data.types.plot.objects$plot,
					num.genes.heatmap = num.genes.plot.objects$plot,
					num.sig.genes.heatmap = num.sig.genes.plot.objects$plot
					),
				legends = list(
					data.types.legend = data.types.plot.objects$legend,
					num.genes.legend = num.genes.plot.objects$legend,
					num.sig.genes.legend = num.sig.genes.plot.objects$legend
					)
			));
		}
	}
