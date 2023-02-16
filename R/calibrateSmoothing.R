#' Automatic calibration of the smoothing parameter
#'
#' @param x A 'mareyMap' object
#' @param method An interpolation method, either 'loess' or 'splines'
#' @param nResampling Number of iterations in the cross-validation procedure
#' @param K Number of clusters to subset in K-fold cross-validation
#' @param from Minimum of the range of smoothing values to explore
#' @param to Maximum of the range of smoothing values to explore
#' @param nCores Number of cores to parallelize
#' @param degree (optional) The degree parameter of the polynomial used in the 'loess' method
#'
#' @return The best smoothing parameter value
#' @export
#'
#' @import pbmcapply
#'
#' @examples calibrateSmoothing(x, "loess")
calibrateSmoothing = function(x, method = "loess",
                              nResampling = 1000, K = 5,
                              from = 0.2, to = 0.5,
                              nCores = 1, degree = 2) {
  if (nrow(x) > 10000) {
    warning("Dataset too large (>10,000 markers). Subsampling 10%...\n")
    x = x[sample(1:nrow(x), size = nrow(x)/10, replace = FALSE),]
  }

  smooth2test = seq(from = from, to = to, by = 0.05)
  if (nrow(x) < 90) {
    warning("Number of markers is low, the space parameter of smoothing is reduced.\n")
    smooth2test = seq(from = max(0.3, from), to = to, by = 0.05)
  }

  size = length(smooth2test)*nrow(x)
  crossvalidation = data.frame(smoothParam = smooth2test, MSE = NA)

  for (s in smooth2test) {
    MSEBoot = numeric(nResampling)
    cat("Fitting smoothing parameter", s, "...\n")
    MSEBoot = unlist(pbmcapply::pbmclapply(X = 1:nResampling,
                                function(X) {cvResampling(x, K = K, smooth = s, method = method, degree = degree)},
                                mc.cores = nCores))
    crossvalidation$MSE[which(crossvalidation$smoothParam == s)] = mean(MSEBoot, na.rm = TRUE)
  }

  smoothParam = crossvalidation$smoothParam[crossvalidation$MSE == min(crossvalidation$MSE, na.rm = TRUE)]

  return(smoothParam)
}




