#' Variance of recombination rates
#'
#' @param x a `marey_map` object.
#'
#' @return (numeric) the variance of recombination rates.
#'
#' @export
variance = function(x) {
  rec = x$recMap$recRate
  m = var(rec, na.rm = TRUE)
  return(m)
}
