#' Fit a smoothing spline function to the Marey map
#'
#' @description
#' Fit a smoothing spline to the Marey map with the 'npreg' statistical package
#' which offers more smoothing parameter selection methods and more smoothing types than the regular 'smooth.spline()'.
#' NOTE that for now the automatic calibration only tune the 'spar' parameter, with other parameters set to DEFAULT.
#'
#' @param x a 'mareyMap' data frame
#' @param spar a smoothing parameter, typically (but not necessarily) in (0,1].
#' @param m Penalty order (integer). The penalty functional is the integrated squared m-th derivative of the function. Defaults to m=2, which is a cubic smoothing spline. Set m=1 for a linear smoothing spline or m=3 for a quintic smoothing spline.
#' @param nknots Positive integer or function specifying the number of knots. Ignored if either all.knots = TRUE or the knot values are input using the knots argument.
#' @param all.knots If TRUE, all distinct points in x are used as knots. If FALSE (default), a sequence knots is placed at the quantiles of the unique x values; in this case, the input nknots specifies the number of knots in the sequence. Ignored if the knot values are input using the knots argument.
#' @param lambda Computational smoothing parameter. This value is weighted by n to form the penalty coefficient (see Details). Ignored if spar is provided.
#'
#' @return a smoothing spline object.
#' @export
#'
#' @import npreg
fit_spline = function(x, spar = numeric(), m = 2, nknots = .nknots.smspl, all.knots = FALSE, lambda = NULL) {
  fit = ss(x = x$phys, y = x$gen, spar = spar, m = m, nknots = nknots, all.knots = all.knots, lambda = lambda)
  return(fit)
}
