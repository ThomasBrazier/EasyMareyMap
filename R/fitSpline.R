#' Fit a cubic smoothing spline function to the Marey map
#'
#' @param x a 'mareyMap' object.
#' @param spar a smoothing parameter, typically (but not necessarily) in (0,1].
#'
#' @return a cubic smoothing spline object.
#' @export
#'
#' @importFrom stats smooth.spline
#'
fitSpline = function(x, spar = numeric()) {
  fit = smooth.spline(x = x$phys, y = x$gen, spar = spar)
  return(fit)
}
