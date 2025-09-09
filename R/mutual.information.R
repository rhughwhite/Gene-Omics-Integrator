#' Define number of bins for data discretization
#'
#' Calculates the number of bins for discretization according to the Freedman–Diaconis rule.
#'
#' @param values continuous data to discretize
#' @return number of bins estimated through the Freedman–Diaconis rule
#' @export
bins.fd <- function(values) {
	values.IQR <- IQR(values);
	if (values.IQR == 0) { # if IQR is 0, n.bins cannot be calculated
		return(NA);
		} else {
		return(diff(range(values)) / (2 * IQR(values) / length(values)^(1 / 3)));
		}
	}

#' Calculate mutual information between two vectors of data
#'
#' Perform discretization and calculation of mutual information using functions from entropy library
#'
#' @param data.1 first input vector from which to calculate mutual information
#' @param data.2 second input vector from which to calculate mutual information
#' @param data.1.categorical does data.1 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values of data.1.
#' data.1 must be numeric for discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param data.2.categorical does data.1 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values of data.1.
#' data.1 must be numeric for discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param n.bins.1 number of bins to use in discretizing data.1. Used if FD.n.bins = FALSE
#' @param n.bins.2 number of bins to use in discretizing data.2. Used if FD.n.bins = FALSE. If not specified, defaults to n.bins.1.
#' @param FD.n.bins whether to calculate number of bins for discretization based upon input data using the Freedman–Diaconis rule.
#' Performed for both input vectors separately. Overrides n.bins and sample.size.bins.
#' @param equal.frequency.bins whether to use an adaptive bin width which aims to equalize the frequencies of observations in each bin. Defaults to FALSE
#' @param normalize.mutual.info should the mutual information be normalized (bounded between 0-1). Applies the 'Symmetric Uncertainty' method of Witten & Frank 2005. Defaults to TRUE.
#' @param use.shrinkage whether to perform shrinkage when calculating mutual information using the James-Stein estimator. Defaults to FALSE
#' @param simple.sample.size.correction apply a naive correction to reduce the bias coming from a limited sample size
#' @return mutual information (after optional normalization)
#' @export
#' @references
#' I. H. Witten, and E. Frank, Data Ming: Practical Machine Learning Tools and Techniques. Morgan Kaufmann, 2nd Edition, 2005.
#'

mutual.information <- function(
	data.1,
	data.2,
	data.1.categorical = FALSE,
	data.2.categorical = FALSE,
	n.bins.1,
	n.bins.2 = NULL,
	FD.n.bins = FALSE,
	equal.frequency.bins = FALSE,
	normalize.mutual.info = TRUE,
	use.shrinkage = FALSE,
	simple.sample.size.correction = FALSE) {

	if (FD.n.bins) {
		n.bins.1 <- bins.fd(data.1);
		n.bins.2 <- bins.fd(data.2);
		if (is.na(n.bins.1) | is.na(n.bins.2)) {
			return(NA); # if IQR is 0, n.bins cannot be calculated
			}
		} else if (is.null(n.bins.2)) {
		n.bins.2 <- n.bins.1;
		}
	if (data.1.categorical) {
		n.bins.1 <- length(unique(data.1));
		}
	if (data.2.categorical) {
		n.bins.2 <- length(unique(data.2));
		}
	# equal frequency bins used only if data are not categorical
	if (equal.frequency.bins & (!data.1.categorical | !data.2.categorical)) {
		if (data.1.categorical) {
			data.binned <- table(
				cut(
					x = data.1,
					breaks = n.bins.1
					),
				Hmisc::cut2(
					x = data.2,
					g = n.bins.2,
					m = length(data.2) / n.bins.2
					)
				);
			} else if (data.2.categorical) {
				data.binned <- table(
					Hmisc::cut2(
						data.1,
						g = n.bins.1,
						m = length(data.1) / n.bins.1
						),
					cut(
						data.2,
						breaks = n.bins.2
						)
					);
			} else {
				data.binned <- table(
					Hmisc::cut2(
						data.1,
						g = n.bins.1,
						m = length(data.1) / n.bins.1
						),
					Hmisc::cut2(
						data.2,
						g = n.bins.2,
						m = length(data.2) / n.bins.2
						)
					);
			}
		} else { #equal width bins
		data.binned <- entropy::discretize2d(data.1, data.2, n.bins.1, n.bins.2);
		}

	if (use.shrinkage) {
		H1 <- entropy::entropy.shrink(rowSums(data.binned), unit = 'log2', verbose = FALSE)[[1]]; #H of marginal probabilities
		H2 <- entropy::entropy.shrink(colSums(data.binned), unit = 'log2', verbose = FALSE)[[1]];
		I <- entropy::mi.shrink(data.binned, unit = 'log2', verbose = FALSE)[[1]];
		} else {
		H12 <- entropy::entropy(data.binned, unit = 'log2'); #H of joint probability
		H1 <- entropy::entropy(rowSums(data.binned), unit = 'log2'); #H of marginal probabilities
		H2 <- entropy::entropy(colSums(data.binned), unit = 'log2');
		I <- H1 + H2 - H12; #calculate I (mutual information)
		}
	if (simple.sample.size.correction) {
		I <- I - ((mean(c(n.bins.1, n.bins.2)) ^ 2) / (2 * length(data.1)));
		}
	if (normalize.mutual.info) {
		I <- (I * 2) / (H2 + H1);
		}
	return(I);
	}

#' Wrapper function - Calculate mutual information
#'
#' Calls mutual.information function. Handles cases with invariate data. Handles optional permutation-based significance estimation
#'
#' @param data.1 first input vector from which to calculate mutual information
#' @param data.2 second input vector from which to calculate mutual information
#' @param data.1.categorical does data.1 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values of data.1.
#' data.1 must be numeric for discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param data.2.categorical does data.1 represent a categorical variable. Used if FD.n.bins is TRUE. If so, the number of bins used is equal to the number of unique values of data.1.
#' data.1 must be numeric for discretization, therefore categorical data must still be encoded in numeric form. Defaults to FALSE.
#' @param n.bins.1 number of bins to use in discretizing data.1. Used if FD.n.bins = FALSE
#' @param n.bins.2 number of bins to use in discretizing data.2. Used if FD.n.bins = FALSE. If not specified, defaults to n.bins.1.
#' @param FD.n.bins whether to calculate number of bins for discretization based upon input data using the Freedman–Diaconis rule.
#' Performed for both input vectors separately. Overides n.bins. Performed for both input vectors separately.
#' If the IQR for a given set of values is 0, the number of bins cannot be calculated and an NA value for mutual information will be returned.
#' @param equal.frequency.bins whether to use an adaptive bin width which aims to equalize the frequencies of observations in each bin. Defaults to FALSE
#' @param normalize.mutual.info should the mutual information be normalized (bounded between 0-1). Defaults to TRUE.
#' @param calculate.pvalue whether to calculate probability of mutual information using a permutation-based approach
#' @param n.permutations if calculate.pvalue = TRUE, determines number of permutations to use. Defaults to 10^5.
#' @param use.shrinkage whether to perform shrinkage when calculating mutual information using the James-Stein estimator. Defaults to FALSE
#' @param simple.sample.size.correction apply a naive correction to reduce the bias coming from a limited sample size
#' @return If calculate.pvalue = FALSE, mutual information is returned (after optional normalization).
#' If calculate.pvalue = TRUE, a vector of length three is returned continaing:
#' observed/true mutual information, expected mutual information = median of purmutations, p.value = probability of observed mutual information
#' @export
get.mutual.information <- function(
	data.1,
	data.2,
	data.1.categorical = FALSE,
	data.2.categorical = FALSE,
	n.bins.1,
	n.bins.2 = NULL,
	FD.n.bins = FALSE,
	equal.frequency.bins = FALSE,
	calculate.pvalue,
	n.permutations = 10^5,
	normalize.mutual.info = TRUE,
	use.shrinkage = FALSE,
	simple.sample.size.correction = FALSE
	) {
	#handle cases of invariate data
	if (var(data.1) == 0 | var(data.2) == 0) {
		I <- 0;
		if (calculate.pvalue) {
			p.value <- NA;
			return(c(I, I, p.value));
		} else {
			return(I);
			}
	} else {
		if (is.null(n.bins.2)) {
			n.bins.2 <- n.bins.1;
			}
		#calculate mutual.information
		I <- mutual.information(
			data.1 = data.1,
			data.2 = data.2,
			data.1.categorical = data.1.categorical,
			data.2.categorical = data.2.categorical,
			n.bins.1 = n.bins.1,
			n.bins.2 = n.bins.2,
			FD.n.bins = FD.n.bins,
			equal.frequency.bins = equal.frequency.bins,
			normalize.mutual.info = normalize.mutual.info,
			use.shrinkage = use.shrinkage,
			simple.sample.size.correction = simple.sample.size.correction
			);
		if (is.na(I)) {# if IQR is 0, n.bins cannot be calculated and a value of NA will be returned for mutual information
			if (calculate.pvalue) {
				p.value <- NA;
				return(c(I, I, p.value));
				} else {
					return(I);
				}
			}
		## Calculating p-value
		if (calculate.pvalue) {

			#calculate randomized/expected mutual.information
			I.exp <- sapply(seq(n.permutations), function(permutation) #enumerate n.permutations (permutation = dummy variable)
				mutual.information(
					data.1 = sample(data.1),
					data.2 = data.2,
					data.1.categorical = data.1.categorical,
					data.2.categorical = data.2.categorical,
					n.bins.1 = n.bins.1,
					n.bins.2 = n.bins.2,
					FD.n.bins = FD.n.bins,
					equal.frequency.bins = equal.frequency.bins,
					normalize.mutual.info = normalize.mutual.info,
					use.shrinkage = use.shrinkage,
					simple.sample.size.correction = simple.sample.size.correction
					)
				);
			#calculate p.value
			p.value <- sum(I.exp >= I) / n.permutations;
			return(c(I, median(I.exp), p.value));
		} else {
			return(I);
			}
		}
	}
