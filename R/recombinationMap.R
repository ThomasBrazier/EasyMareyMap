#' Estimate the recombination map
#'
#' @param x A 'mareyMap' object
#' @param chromosome The name of the chromosome to process
#' @param method An interpolation method, either 'loess' or 'splines'
#' @param K Number of clusters to subset in K-fold cross-validation
#' @param boot Number of bootstrap to estimate a 95% confidence interval
#' @param nCores Number of cores to parallelize
#' @param smoothing (optional) Smoothing parameter, if you do not want to calibrate it
#' @param nResampling Number of iterations in the cross-validation procedure
#' @param calibrationRange Range of the smoothing parameter space to explore
#' @param degree (optional) The degree parameter of the polynomial used in the 'loess' method
#' @param windowSize Size of the windows along the chromosome in basepairs (default = 10^5)
#'
#' @return A 'mareyMap' object with the estimated recombination map
#' @export
#'
#' @import pbmcapply
#' @import parallel
#'
#' @examples recombinationMap(x, chromosome = "1", method = "loess")
recombinationMap = function(x, chromosome = character(), method = c("loess", "spline"),
                            K = 5, boot = 1000, nCores = 1,
                            smoothing = numeric(), nResampling = 1000,
                            calibrationRange = c(0.2, 0.5),
                            degree = 2, windowSize = 10^5) {
  if (!(method %in% c("loess", "spline"))) {
    stop("Interpolation method must be specified: loess or spline.")
  }

  calibrationFrom = calibrationRange[1]
  calibrationTo = calibrationRange[2]

  chromosome = as.character(chromosome)

  nCores = min(nCores,
               (parallel::detectCores(all.tests = FALSE, logical = TRUE)-1))

  cat("Processing chromosome ", chromosome, "...\n")
  df = x$mareyMap
  df = df[which(df$map == chromosome), ]
  df = df[!is.na(df$gen), ]
  df = df[!is.na(df$phys), ]

  if (length(smoothing) > 0) {
    if (is.numeric(smoothing)) {
      smoothingParam = smoothing
    }
  } else {
    cat("Calibrate the smoothing parameter.\n")
    smoothingParam = calibrateSmoothing(df,
                                        method = method,
                                        nResampling = nResampling,
                                        K = K,
                                        from = calibrationFrom,
                                        to = calibrationTo,
                                        nCores = nCores,
                                        degree = degree
    )
  }

  x$interpolationMethod[chromosome] = method
  x$smoothingParam[chromosome] = smoothingParam
  x$windowSize[chromosome] = windowSize
  x$nBootstrap[chromosome] = boot

  cat("Fit the model.\n")
  boundaries = seq(0, max(df$phys, na.rm = TRUE), by = windowSize)
  physWindows = data.frame(start = boundaries[-length(boundaries)],
                           end = (boundaries[-1] - 1)) # A sequence of windows coordinates
  physWindows$point = round((physWindows$end + physWindows$start)/2, digits = 0)
  bootPrediction = matrix(NA, nrow = nrow(physWindows) - 1, ncol = boot)

  pb = txtProgressBar(min = 1, max = boot, initial = 1)
  for (b in 1:boot) {
    setTxtProgressBar(pb, b)

    resampled = df[sample(1:nrow(df), replace = TRUE), ]
    if (method == "loess") {
      fitMarey = fitLoess(resampled, span = smoothingParam, degree = degree)
      predicted = predict(fitMarey, newdata = physWindows$point)
    }
    if (method == "spline") {
      fitMarey = fitSpline(resampled, spar = smoothingParam)
      predicted = predict(fitMarey, physWindows$point)$y
    }
    bootPrediction[,b] = (diff(predicted) / diff(physWindows$point))
  }
  Sys.sleep(1)
  close(pb)

  cat("Estimate recombination rates.\n")
  dX = rowMeans(embed(physWindows$point, 2))
  upperCI = function(x) {
    quantile(x, 0.975, na.rm = TRUE)
  }
  lowerCI = function(x) {
    quantile(x, 0.025, na.rm = TRUE)
  }
  dY = rowMeans(bootPrediction, na.rm = TRUE)
  dYupper = apply(bootPrediction, MARGIN = 1, FUN = upperCI)
  dYlower = apply(bootPrediction, MARGIN = 1, FUN = lowerCI)
  dY[dY < 0] = 0
  dYupper[dYupper < 0] = 0
  dYlower[dYlower < 0] = 0

  estimates = data.frame(
    set = unique(df$set),
    map = chromosome,
    start = dX - windowSize,
    end = dX + windowSize,
    recRate = dY,
    upperRecRate = dYupper,
    lowerRecRate = dYlower
  )

  x$recMap = x$recMap[x$recMap$map != chromosome,]
  x$recMap = rbind(x$recMap, estimates)

  return(x)
}