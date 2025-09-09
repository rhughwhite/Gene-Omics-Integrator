#' Calculate mutual information between rows of two matrices
#'
#' Helper function - runs get.MI on two matrices with common rows (genes)
#' Matching of rows and colnames is optionally performed
#' Matching colnames (i.e. samples) and row.names (i.e. genes) are identified and used to order the two input matrices
#' @param input.matrix.1 first input matrix from which to calculate mutual information
#' @param input.matrix.2 second input matrix from which to calculate mutual information
#' @param input.matrix.1.categorical does input.matrix.1 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values for each row of input.matrix.1 separately.
#' input.matrix.1 must be numeric to allow downstream discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param input.matrix.1.categorical does input.matrix.2 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values for each row of input.matrix.2 separately.
#' input.matrix.2 must be numeric to allow downstream discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param match.rows if TRUE, common row.names of input.matrices are identified and used for computation of mutual information. Defaults to TRUE.
#' If FALSE, input matrices should have the same number of rows, which will be compared in order i.e. the first row of input.matrix.1 compared with the first row of input.matrix.2. Defaults to TRUE.
#' @param match.columns if TRUE, columns of input matrices will be orderd by shared colnames (i.e. sample names). If FALSE, columns of input matrices should already be matched.
#' @param filter.rows whether to filter rows from each data type pair after intersecting by column names (samples). This may be desirable as data may become 'too sparse' after intersecting
#' between a given data type pair. Defaults to FALSE.
#' @param filter.threshold.1 used if filter.rows is TRUE. What is the minimum proportion of samples that should have non-zero values for each row of input.matrix.1.
#' e.g. 0.5 would indicate >= 50% of samples must have non-zero observations.
#' @param filter.threshold.2 used if filter.rows is TRUE. What is the minimum proportion of samples that should have non-zero values for each row of input.matrix.2.
#' @param normalize.input should data be normalized using the inverse normal transformation.
#' See `?RNOmni::RankNorm`. Normalization is performed after intersecting data types by sample and optional filtering of rows. This may be more desirable than normalizing data before intersection and sampling.
#' Defaults to FALSE
#' @param n.bins number of bins to use in discretizing data, applied to both input vectors. Used if FD.n.bins = FALSE
#' @param FD.n.bins whether to calculate number of bins for discretization based upon input data using the Freedman–Diaconis rule.
#' Performed for both input vectors separately. Overides n.bins and sample.size.bins. Performed for both input vectors separately. Defaults to FALSE
#' @param sample.size.bins whether to calculate number of bins for discretization as a function of sample size (N).
#' Bin number = N^(1/3). Defaults to FALSE
#' @param equal.frequency.bins whether to use an adaptive bin width which aims to equalize the frequencies of observations in each bin. Defaults to FALSE
#' @param normalize.mutual.info should the mutual information be normalized (bounded between 0-1)
#' @param calculate.pvalue whether to calculate probability of mutual information using a permutation-based approach. Defaults to FALSE
#' @param n.cores how many cores to use for parallelization of mutual information calculation. Parallelization occurs across genes/rows of input matrices. Defaults to 1.
#' @param n.permutations if calculate.pvalue = TRUE, determines number of permutations to use
#' @param use.shrinkage whether to perform shrinkage when calculating mutual information using the James-Stein estimator. Defaults to FALSE
#' @param simple.sample.size.correction apply a naive correction to reduce the bias coming from a limited sample size
#' @return matrix containing mutual information for every matched row of input matrices.
#' If calculate.pvalue = TRUE, three columns are returned:
#' observed/true mutual information, expected mutual information = median of permutations, p.value = probability of observed mutual information
#' the number of rows will be the length of the intersections of the row.names between input matrices
#' @importFrom foreach %dopar%
#' @export
get.mutual.information.from.matrices <- function(
	input.matrix.1,
	input.matrix.2,
	input.matrix.1.categorical = FALSE,
	input.matrix.2.categorical = FALSE,
	match.rows = TRUE,
	match.columns = TRUE,
	filter.rows = FALSE,
	filter.threshold.1,
	filter.threshold.2,
	normalize.input = FALSE,
	normalize.mutual.info = TRUE,
	n.permutations,
	calculate.pvalue = FALSE,
	n.cores = 1,
	n.bins,
	FD.n.bins = FALSE,
	sample.size.bins,
	equal.frequency.bins,
	use.shrinkage = FALSE,
	simple.sample.size.correction = FALSE
	) {

	if (match.columns) {
		#extract the samples which are common to data type pair
		samples <- intersect(colnames(input.matrix.1), colnames(input.matrix.2));
		input.matrix.1 <- input.matrix.1[, samples];
		input.matrix.2 <- input.matrix.2[, samples];
		}

	if (filter.rows) {
		input.matrix.1 <- input.matrix.1[which(apply(
			input.matrix.1,
			1,
			function(this.row)
				length(which(this.row != 0)) >= round(filter.threshold.1 * ncol(input.matrix.1))
			)), ];
		input.matrix.2 <- input.matrix.2[which(apply(
			input.matrix.2,
			1,
			function(this.row)
				length(which(this.row != 0)) >= round(filter.threshold.2 * ncol(input.matrix.2))
			)), ];
		}

	if (match.rows) {
		#extract the rows which are common to data type pair
		common.rows <- intersect(row.names(input.matrix.1), row.names(input.matrix.2));
		input.matrix.1 <- input.matrix.1[common.rows, ];
		input.matrix.2 <- input.matrix.2[common.rows, ];
		}
	#if no data left return dummy results
	if (0 == nrow(input.matrix.1) | 0 == nrow(input.matrix.1)) {
		return(NULL);
		}

	if (normalize.input) {
		input.matrix.1 <- normalize.data(input.matrix.1);
		input.matrix.2 <- normalize.data(input.matrix.2);
		}
	#calculate n.bins based on sample size if user specified
	if (sample.size.bins) {
		n.bins <- floor(ncol(input.matrix.1) ^ (1 / 3));
		}
	#calculate mutual information results
	cluster <- parallel::makeCluster(n.cores);
	doParallel::registerDoParallel(cores = cluster);
	MI.results <- foreach::foreach(this.row = seq(nrow(input.matrix.1)), .combine = 'rbind') %dopar% {
		get.mutual.information(
			data.1 = input.matrix.1[this.row, ],
			data.2 = input.matrix.2[this.row, ],
			data.1.categorical = input.matrix.1.categorical,
			data.2.categorical = input.matrix.2.categorical,
			n.permutations = n.permutations,
			calculate.pvalue = calculate.pvalue,
			FD.n.bins = FD.n.bins,
			equal.frequency.bins = equal.frequency.bins,
			n.bins.1 = n.bins,
			normalize.mutual.info = normalize.mutual.info,
			use.shrinkage = use.shrinkage,
			simple.sample.size.correction = simple.sample.size.correction
			);
		}
	parallel::stopCluster(cluster);
	#add dimnames to results
	if (calculate.pvalue) {
		colnames(MI.results) <- c('MI', 'expected.MI', 'p.value');
		} else {
		MI.results <- matrix(MI.results, dimnames = list(NULL, 'MI'));
		}
	if (match.rows) {
		row.names(MI.results) <- common.rows;
		}
	return(MI.results);

	on.exit({#TODO - test this
		try({
			cat('Attempting to stop cluster\n');
			doParallel::stopImplicitCluster(cluster);
			})
		})
	}

#' Compute mutual information from list of data.frames/matrices
#'
#' Rows of data from each data.frame/matrix will be used to compute mutual information. Each pairwise combination of inputs will be compared, excluding self-comparisons
#' Either rows will be matched for computation of mutual information via common row.names (if match.rows is TRUE), or each row pair of rows will be matched in turn
#' Matching of rows and colnames is optionally performed
#' Matching colnames (i.e. samples) and row.names (i.e. genes) are identified between pairs of input data.frames/matrices
#' @param omics.data.list a named list of data.frames/matrices with samples in columns and gene-wise data in rows. Names will be used for detailing pariwse comparisons in output.
#' @param data.types.categorical a named vector or list specifying whether each data type in omics.data represents categorical data. Used only if FD.n.bins is TRUE.
#' If so, the number of bins used is equal to the number of unique values for each row of input.matrix.1 separately. Defaults to FALSE for all elements of omics.data.
#' @param comparisons a named list that details the pairwise comparisons to be made. Each element of the list should be a vector of length two indicating a pair of data types from omics.data.list which will be compared.
#' Names of the list will be used for outputs. Defaults to NULL (compare all input data in a pairwise fashion)
#' @param match.rows if TRUE, common row.names of input.matrices are identified and used for computation of mutual information
#' If FALSE, input matrices should have the same number of rows, which will be compared in order i.e. the first row of input.matrix.1 compared with the first row of input.matrix.2. Defaults to TRUE.
#' @param match.columns if TRUE, columns of input matrices will be orderd by shared colnames (i.e. sample names). If FALSE, columns of input matrices should already be matched. Defaults to TRUE.
#' @param filter.rows whether to filter rows from each data type pair after intersecting by column names (samples). This may be desirable as data may become 'too sparse' after intersecting
#' between a given data type pair. Defaults to FALSE.
#' @param filter.threshold used if filter.rows is TRUE. What is the minimum proportion of samples that should have non-zero values for each data type pair after intersecting by sample. e.g. 0.5 would indicate >= 50% of samples must have non-zero observations.
#' This can be a single value for all data types, or a named vector/list with a threshold for each data type.
#' @param normalize.input should data be normalized using the inverse normal transformation.
#' See `?RNOmni::RankNorm`. Normalization is performed after intersecting data types by sample and optional filtering of rows. This may be more desirable than normalizing data before intersection and sampling.
#' Defaults to FALSE
#' @param n.bins number of bins to use in discretizing data, applied to both input vectors. Used if FD.n.bins = FALSE.
#' @param FD.n.bins whether to calculate number of bins for discretization based upon input data using the Freedman–Diaconis rule.
#' Performed for both input vectors separately. Overides n.bins and sample.size.bins. Performed for both input vectors separately. Defaults to FALSE
#' @param sample.size.bins whether to calculate number of bins for discretization as afunction of sample size (N).
#' Bin number = N^(1/3). Defaults to FALSE
#' @param equal.frequency.bins whether to use an adaptive bin width which aims to equalize the frequencies of observations in each bin. Defaults to FALSE
#' @param normalize.mutual.info should the mutual information be normalized (bounded between 0-1)
#' @param calculate.pvalue whether to calculate probability of mutual information using a permutation-based approach. Defaults to FALSE
#' @param FDR.threshold threshold for determining number of significant genes in returned summary statistics data.frame. Defaults to 0.05
#' @param n.cores how many cores to use for parallelization of mutual information calculation. Parallelization occurs across genes/rows of input matrices. Defaults to 1.
#' @param n.permutations if calculate.pvalue = TRUE, determines number of permutations to use
#' @param use.shrinkage whether to perform shrinkage when calculating mutual information using the James-Stein estimator. Defaults to FALSE
#' @param simple.sample.size.correction apply a naive correction to reduce the bias coming from a limited sample size
#' @return A named list with an entry for every pairwise comparison of input matrices/data.frames (excluding self comparisons)
#' If calculate.pvalue = FALSE, a single matrix containing mutual information is returned
#' If calculate.pvalue = TRUE, a list of 3 matrices is returned with: mutual information, expected (median of permuted) mutual information, p value probability of observed mutual information
#' @export
#' @examples
#' data('simple.data');
#'
#'data.types.categorical <- c(
#'    data.a = FALSE,
#'    data.b = FALSE,
#'    data.c = TRUE # data.c contains categorical data
#'    );
#'
#'### calculate mutual information between all pairs of datasets
#'set.seed(123);
#'mutual.info <- calculate.mutual.information.pairwise.data.types(
#'    omics.data.list = simple.data,
#'    data.types.categorical = data.types.categorical,
#'    match.rows = TRUE,
#'    match.columns = TRUE,
#'    sample.size.bins = TRUE,
#'    calculate.pvalue = TRUE,
#'    n.permutations = 100
#'    );
#'
#'### plot results
#'data.type.colours <- setNames(BoutrosLab.plotting.general::default.colours(length(simple.data)), names(simple.data));
#'
#'data.type.legend.labels <- names(data.type.colours);
#'
#'miplot <- create.mutual.information.summary.plot(
#'    MI = mutual.info$MI,
#'    order.by.median = TRUE,
#'    data.type.colours = data.type.colours,
#'    data.type.legend.labels = data.type.legend.labels,
#'    summary.statistics = mutual.info$summary.statistics,
#'    return.plot = TRUE, # to save to file, set FALSE and use filename = '<plot-name>.tiff'
#'    return.plot.elements = FALSE,
#'    plot.width = 8,
#'    plot.height = 5
#'    );
#'miplot;
calculate.mutual.information.pairwise.data.types <- function(
	omics.data.list,
	data.types.categorical = setNames(rep(FALSE, length(omics.data.list)), names(omics.data.list)),
	comparisons = NULL,
	match.rows = TRUE,
	match.columns = TRUE,
	filter.rows = FALSE,
	filter.threshold = NULL,
	normalize.input = FALSE,
	n.permutations,
	calculate.pvalue = FALSE,
	FDR.threshold = 0.05,
	n.cores = 1,
	n.bins = NULL,
	FD.n.bins = FALSE,
	sample.size.bins = TRUE,
	equal.frequency.bins = FALSE,
	normalize.mutual.info = TRUE,
	use.shrinkage = FALSE,
	simple.sample.size.correction = FALSE
	) {

	#cast to matrix for more efficient manipulation
	omics.data.list <- lapply(omics.data.list, function(data.type) as.matrix(data.type));

	if (is.null(comparisons)) { #make a named list that details the pairwise comparisons to be made
		comparisons <- combn(names(omics.data.list), 2);
		colnames(comparisons) <- apply(comparisons, 2, paste, collapse = '_');
		comparisons <- unlist(apply(comparisons, 2, list), recursive = FALSE);
		}
	#prepare filter.threshold if required
	if (filter.rows & 1 == length(filter.threshold)) {
		filter.threshold <- setNames(rep(filter.threshold, length(omics.data.list)), names(omics.data.list));
		}
	#run the mutual.information functions for each of the pairwise comparisons
	MI.results <- lapply(comparisons,
		function(compare)
			get.mutual.information.from.matrices(
				input.matrix.1 = omics.data.list[[compare[1]]],
				input.matrix.2 = omics.data.list[[compare[2]]],
				input.matrix.1.categorical = data.types.categorical[[compare[1]]],
				input.matrix.2.categorical = data.types.categorical[[compare[2]]],
				match.rows = match.rows,
				match.columns = match.columns,
				filter.rows = filter.rows,
				filter.threshold.1 = filter.threshold[[compare[1]]],
				filter.threshold.2 = filter.threshold[[compare[2]]],
				n.permutations = n.permutations,
				calculate.pvalue = calculate.pvalue,
				n.cores = n.cores,
				FD.n.bins = FD.n.bins,
				sample.size.bins = sample.size.bins,
				equal.frequency.bins = equal.frequency.bins,
				normalize.mutual.info = normalize.mutual.info,
				n.bins = n.bins,
				use.shrinkage = use.shrinkage,
				simple.sample.size.correction = simple.sample.size.correction
				)
		);
	#filter missing data - present if no row/gene intersection between a given data pair
	MI.results <- Filter(function(data.pair) !is.null(data.pair), MI.results);
	#add specific names to each pairwise comparison
	MI.results <- lapply(names(MI.results),
		function(compare)
			setNames(
				data.frame(MI.results[[compare]]),
				paste(compare, colnames(MI.results[[compare]]), sep = '_')
				)
		);
	#merge pairwise comparison
	MI.results <- Reduce(function(x, y)
		transform(
			merge(x, y, by = 0, all = TRUE),
			row.names = Row.names, Row.names = NULL
			),
		MI.results
		);
	#separate out each data type and return i.e. turn 1 matrix into 3 seperate matrices
	if (calculate.pvalue) {
		results <- list(
			MI = setNames(
				MI.results[, grep('_MI$', colnames(MI.results)), drop = FALSE],
				sub('_MI', '', grep('_MI$', colnames(MI.results), value = TRUE))
				),
			expected.MI = setNames(
				MI.results[, grep('expected.MI$', colnames(MI.results)), drop = FALSE],
				sub('_expected.MI', '', grep('_expected.MI$', colnames(MI.results), value = TRUE))
				),
			p.value = setNames(
				MI.results[, grep('p.value$', colnames(MI.results)), drop = FALSE],
				sub('_p.value', '', grep('_p.value$', colnames(MI.results), value = TRUE))
				)
			);
		#add back missing data - present if no row/gene intersection between a given data pair
		results <- lapply(results, function(result) {
			result[names(comparisons)[!names(comparisons) %in% names(result)]] <- NA;
			return(result);
			});
		results <- lapply(results, function(result)
			result[, names(comparisons), drop = FALSE]
			);

		#calculate FDR
		results$FDR <- apply(results$p.value, 2, p.adjust, method = 'fdr');

		#generate summary statistics
		num.genes <- apply(results$p.value, 2, function(data.type) length(which(!is.na(data.type))));
		num.sig.genes <- apply(results$FDR, 2, function(data.type) length(which(data.type < FDR.threshold)));
		median.MI <- apply(results$MI, 2, median, na.rm = TRUE);
		expected.median.MI <- apply(results$expected.MI, 2, median, na.rm = TRUE);
		summary.statistics <- data.frame(
			comparison = names(median.MI),
			median.MI = median.MI,
			expected.median.MI = expected.median.MI,
			num.genes = num.genes,
			num.sig.genes = num.sig.genes
			);
		results$summary.statistics <- summary.statistics;
		return(results);
	} else {
		results <- list(MI = setNames(MI.results, sub('_MI', '', colnames(MI.results))));
		results$MI[names(comparisons)[!names(comparisons) %in% names(results$MI)]] <- NA;
		results$MI <- results$MI[, names(comparisons)];
		return(results);
		}

	}

#' Normalize data using the rank-based inverse normal transform.
#'
#' Rows of a matrix of normalized using the rank-based inverse normal transform.
#' See `?RNOmni::RankNorm`
#' @param input.matrix numeric matrix, rows of which will be normalized.
#' @return normalized matrix
#' @export
normalize.data <- function(input.matrix) {
	input.matrix.norm <- t(apply(input.matrix, 1, function(this.row)
		this.row[!is.na(this.row)] <- RNOmni::RankNorm(this.row[!is.na(this.row)])
		));
	return(input.matrix.norm);
	}
