#' @description
#' Interpolate genetic positions (i.e. cumulative gnetic distances)
#' on a new set of physical positions
#'
#' @param phys a vector of physical positions.
#' @param gen a vector of genetic positions.
#' @param new a vector of new physical positions.
#' @param method an interpolation method, either 'loess' or 'splines'.
#' @param smoothing (optional) smoothing parameter, if you do not want to calibrate it.
#' @param degree (optional) the degree parameter of the polynomial used in the 'loess' method.
#'
#' @return a list with 'newGen', a vector of new genetic positions for each new marker position
#' and 'se' the standard error of the interpolated result (if relevant for the method).
#' @export
#'
#' @import pbmcapply
#' @import parallel
#' @import utils
#' @import stats
predictGeneticMap = function(phys,
                             gen,
                             new,
                             method = "loess",
                             smoothing = 0.3,
                             degree = 2) {

  # Fit the Marey function to all points
  if (method == "loess") {
    fitMarey = loess(gen ~ phys, span = smoothing, degree = degree)
  }
  if (method == "spline") {
    fitMarey = smooth.spline(x = phys, y = gen, spar = smoothing)
  }
  
  # Predict new positions fro the fitted model
  if (method == "loess") {
    pred = predict(fitMarey, newdata = new, se = TRUE)
    newGen = pred$fit
    se = pred$se
  }
  if (method == "spline") {
    pred = predict(fitMarey, new)
    newGen = pred$y
    se = NA
  }
  
  newGen[which(newGen < 0)] = 0
  
  return(list(newGen = newGen, se = se))
}
