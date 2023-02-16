#' K-fold systematic cross-validation of the smoothing parameter
#'
#' @param x A marey map (the 'mareyMap' slot of a mareyMap object)
#' @param K Number of clusters to subset in K-fold cross-validation
#' @param smooth Smoothing parameter
#' @param method An interpolation method, either 'loess' or 'splines'
#' @param degree (optional) Degree value of the polynomial function of 'loess' method (default 2)
#'
#' @return Goodness-of-fit criterion (mean squared difference)
#' @export
#'
#' @examples cvResampling(x)
cvResampling = function(x, K = 5, smooth = numeric(), method = "loess", degree = 2) {
  idx = sample(1:nrow(x), size = nrow(x)*((K-1)/K), replace = FALSE)
  training = x[idx,]
  validation = x[-idx,]

  if (method == "loess") {
    fitTraining = fitLoess(training, span = smooth, degree = degree)
    predict = predict(fitTraining, newdata = validation$phys, se = TRUE)
    validation$pred = predict$fit
  }
  if (method == "spline") {
    fitTraining = fitSpline(training, spar = smooth)
    predict = predict(fitTraining, validation$phys)
    validation$pred = predict$y
  }

  df = validation[,c(4,5,7)]
  colnames(df) = c("phys", "obs", "pred")
  mean_diff = mean((df$obs - df$pred)^2, na.rm = TRUE)
  return(mean_diff)
}
