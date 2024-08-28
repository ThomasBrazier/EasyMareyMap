#' Mean of recombination rates
#'
#' @param x a Marey map object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return numeric, the mean of recombination rates.
#'
#' @method mean mareyMap
#' @export
mean.mareyMap = function(x, ...) {
  rec = x$recMap$recRate
  m = mean(rec, na.rm = TRUE)
  return(m)
}
