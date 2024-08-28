#' Median of recombination rates
#'
#' @param x a Marey map object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the median of recombination rates.
#'
#' @method median mareyMap
#' @export
median.mareyMap = function(x, ...) {
  rec = x$recMap$recRate
  m = median(rec, na.rm =  TRUE)
  return(m)
}
