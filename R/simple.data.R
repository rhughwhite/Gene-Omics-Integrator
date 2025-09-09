#' simple.data
#'
#' Simple dataset used for testing package.
#'
#' @format a list of 3 data.frames each with the same features in rows and sample sin columns:
#' \describe{
#' \item{data.a}{continuous features}
#' \item{data.b}{continuous features}
#' \item{data.c}{categorical features}
#' }
#' @examples
#' data('simple.data');
#' lapply(simple.data, dim);
#' lapply(simple.data, function(x) x[1:5,1:5]);
#'
'simple.data'
