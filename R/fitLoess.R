#' Fit a LOESS function to the Marey map
#'
#' @param x a 'mareyMap' data frame.
#' @param span the span parameter which controls the degree of smoothing.
#' @param degree the degree of the polynomials to be used, normally 1 or 2.
#'
#' @return a LOESS object.
#' @export
#'
#' @importFrom stats loess
#'
fitLoess = function(x, span = numeric(), degree = 2) {
  fit = loess(gen ~ phys, x, span = span, degree = degree)
  return(fit)
}
