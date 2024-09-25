#' Interpolate new genetic positions
#'
#' @description
#' Interpolate genetic positions (i.e. cumulative genetic distances)
#' on a new set of physical positions
#' and return a new `marey_map` object with the new set of markers (i.e. new genetic and physical distances)
#' When the `predict_genetic_map` function is used, the standard error of predictions is saved in the 'predict.se' column of the mareyMap slot.
#' (only for the "loess" regression)
#'
#' @param x a `marey_map` object, with a recombination map.
#' @param new a vector of new physical positions.
#'
#' @return a new `marey_map` object with the new set of markers
#' @export
#'
#' @import utils
#' @import stats
#' @importFrom methods is
#' 
predict_genetic_map = function(x, new) {
  # Predict new positions fro the fitted model
  if (is(x$model, "loess")) {
    pred = predict(x$model, newdata = new, se = TRUE)
    newGen = pred$fit
    se = pred$se
  }
  if (is(x$model, "ss")) {
    pred = predict(x$model, new)
    newGen = pred$y
    se = NA
  }
  newGen[which(newGen < 0)] = 0

  newMarey = data.frame(set = unique(x$mareyMap$set),
                        map = unique(x$mareyMap$map),
                        mkr = paste0("predicted_marker_", seq(1, length(new))),
                        phys = new,
                        gen = newGen,
                        vld = TRUE,
                        predict.se = se)
  x$mareyMap = newMarey

  return(x)
}
