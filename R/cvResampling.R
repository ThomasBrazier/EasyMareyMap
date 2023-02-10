cvResampling = function(x, K = 5, smooth = numeric(), method = "loess", degree = 2) {
  idx = sample(1:nrow(x), size = nrow(x)*((K-1)/K), replace = FALSE)
  training = x[idx,]
  validation = x[-idx,]

  if (method == "loess") {
    fitTraining = fitLoess(training, span = s, degree = degree)
    predict = predict(fitTraining, newdata = validation$phys, se = TRUE)
    validation$pred = predict$fit
  }
  if (method == "spline") {
    fitTraining = fitSpline(training, spar = s)
    predict = predict(fitTraining, validation$phys)
    validation$pred = predict$y
  }

  df = validation[,c(4,5,7)]
  colnames(df) = c("phys", "obs", "pred")
  mean_diff = mean((df$obs - df$pred)^2, na.rm = TRUE)
  return(mean_diff)
}
