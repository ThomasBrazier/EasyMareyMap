#' Mean of recombination rates
#'
#' @param x a Marey map object.
#'
#' @return numeric, the mean of recombination rates.
#'
#' @exportS3Method base::mean
#'
#'
mean.mareyMap = function(x) {
  rec = x$recMap$recRate
  m = mean(rec, na.rm = TRUE)
  return(m)
}
