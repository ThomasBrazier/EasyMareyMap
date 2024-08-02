#' Median of recombination rates
#'
#' @param x a Marey map object.
#'
#' @return numeric, the median of recombination rates.
#' @export
#'
#'
median.mareyMap = function(x) {
  rec = x$recMap$recRate
  m = median(rec, na.rm =  TRUE)
  return(m)
}
