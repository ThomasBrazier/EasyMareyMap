#' Weighted Mean of recombination rates
#'
#' @param x a Marey map object.
#'
#' @return numeric, the weighted mean of recombination rates (weighted by windows size).
#' @export
#'
#'
weighted.mean.mareyMap = function(x) {
  rec = x$recMap$recRate
  w = x$recMap$end - x$recMap$start
  m = weighted.mean(rec, w, na.rm = TRUE)
  return(m)
}
