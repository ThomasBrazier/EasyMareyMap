#' Calculate the variance and coefficient of variation of recombination rates
#'
#' @description
#' The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean
#' CV = standard deviation/mean. It is estimated on recombination rates.
#'
#' @param x a `marey_map` object, with a recombination map already estimated
#' 
#' @importFrom stats var
#' @importFrom stats sd
#' 
#' @return a list of numeric values of the variance, standard deviation and coefficient of variation.
#' @export
#'
coefficient_variation = function(x) {
  df = x$recMap

  v = var(df$recRate, na.rm = TRUE)
  std = sd(df$recRate, na.rm = TRUE)
  m = mean(df$recRate, na.rm = TRUE)
  cv = std/m

  return(cv)
}
