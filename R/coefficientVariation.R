#' Calculate the variance and coefficient of variation of recombination rates
#'
#' @description
#' The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean
#' CV = standard deviation/mean. It is estimated on recombination rates.
#'
#' @param x a Marey map object, with a recombination map already estimated.
#'
#' @return a list of numeric values of the variance, standard deviation and coefficient of variation.
#' @export
#'
#'
coefficientVariation = function(x) {
  df = x$recMap

  v = var(df$recRate, na.rm = TRUE)
  std = sd(df$recRate, na.rm = TRUE)
  m = mean(df$recRate, na.rm = TRUE)
  cv = std/m

  return(list(variance = v,
         std.deviation = std,
         CV = cv))
}
