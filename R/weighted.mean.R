#' Weighted Mean of recombination rates
#'
#' @param x a Marey map object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the weighted mean of recombination rates (weighted by windows size).
#'
#' @method weighted.mean mareyMap
#' @export
weighted.mean.mareyMap = function(x, ...) {
  rec = x$recMap$recRate
  w = x$recMap$end - x$recMap$start
  m = weighted.mean(rec, w, na.rm = TRUE)
  return(m)
}
