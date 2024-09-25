#' K-fold systematic cross-validation of the smoothing parameter
#'
#' @param x a `mareyMap` slot of a `marey_map` object).
#' @param K number of clusters to subset in K-fold cross-validation.
#' @param smooth the smoothing parameter.
#' @param method an interpolation method, either 'loess' or 'splines'.
#' @param degree (optional) the degree value of the polynomial function of 'loess' method (default 2).
#'
#' @return the goodness-of-fit criterion (mean squared difference).
#' @export
#'
cv_resampling = function(x, K = 5, smooth = numeric(), method = "loess", degree = 2) {
  idx = sample(1:nrow(x), size = nrow(x)*((K-1)/K), replace = FALSE)
  training = x[idx,]
  validation = x[-idx,]

  if (method == "loess") {
    fitTraining = fit_loess(training, span = smooth, degree = degree)
    predict = predict(fitTraining, newdata = validation$phys, se = TRUE)
    validation$pred = predict$fit
  }
  if (method == "spline") {
    fitTraining = fit_spline(training, spar = smooth)
    predict = predict(fitTraining, validation$phys)
    validation$pred = predict$y
  }

  df = validation[, c("phys", "gen", "pred")]
  colnames(df) = c("phys", "obs", "pred")
  mean_diff = mean((df$obs - df$pred)^2, na.rm = TRUE)
  return(mean_diff)
}
