#' Calculate the Gini index of a recombination map
#'
#' @description
#' The Gini index is a measure of statistical dispersion intended to represent inequality between measures, here recombination rates
#' A Gini coefficient of 0 reflects perfect equality,
#' while a Gini coefficient of 1 reflects maximal inequality among values
#' The Gini index is estimated on recombination rates with the `gini.wtd` function from the `dineq` package.
#'
#' @param x a Marey map object, with a recombination map already estimated.
#' @param bootstrap a numeric value of the number of iterations for a bootstrapped 95\% C.I. Set NULL if no bootstrap.
#'
#' @return a single value of the Gini index if no bootstrap, otherwise a list of numeric values of the Gini index, and lower/upper 95\% C.I.
#' @export
#'
#' @import dineq
gini = function(x, bootstrap = NULL) {
  df = x$recMap

  g = gini.wtd(df$recRate, weights = NULL)

  if (is.null(bootstrap)) {
    return(g)
  } else {

    return(list(gini = g,
                upper = up,
                lower = lo))
  }
}
