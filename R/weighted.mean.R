#' Weighted Mean of recombination rates
#'
#' @param x a `marey_map` object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the weighted mean of recombination rates (weighted by windows size).
#'
#' @importFrom stats weighted.mean
#' @method weighted.mean marey_map
#' @export
weighted.mean.marey_map = function(x, ...) {
  rec = x$recMap$recRate
  w = x$recMap$end - x$recMap$start
  m = weighted.mean(rec, w, na.rm = TRUE)
  return(m)
}
