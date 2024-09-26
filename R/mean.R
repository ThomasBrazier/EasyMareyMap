#' Mean of recombination rates
#'
#' @param x a `marey_map` object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the mean of recombination rates.
#'
#' @method mean marey_map
#' @export
#' 
mean.marey_map = function(x, ...) {
  rec = x$recMap$recRate
  m = mean(rec, na.rm = TRUE)
  return(m)
}
