#' Median of recombination rates
#'
#' @param x a `marey_map` object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the median of recombination rates.
#'
#' @importFrom stats median
#' @method median marey_map
#' @export
#' 
median.marey_map = function(x, ...) {
  rec = x$recMap$recRate
  m = median(rec, na.rm =  TRUE)
  return(m)
}
