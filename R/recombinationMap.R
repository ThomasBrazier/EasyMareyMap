#' Estimate the recombination map
#'
#' @param x a 'mareyMap' object.
#' @param chromosome the name of the chromosome to process.
#' @param method an interpolation method, either 'loess' or 'splines'.
#' @param K number of clusters to subset in K-fold cross-validation.
#' @param boot number of bootstraps to estimate the confidence interval.
#' @param nCores number of cores to parallelize.
#' @param smoothing (optional) smoothing parameter, if you do not want to calibrate it.
#' @param nResampling number of iterations in the cross-validation procedure.
#' @param calibrationRange range of the smoothing parameter space to explore.
#' @param degree (optional) the degree parameter of the polynomial used in the 'loess' method.
#' @param windows size or coordinates of the windows along the chromosome: either integer in basepairs (default = 10^5) or a data.frame of start/end coordinates.
#'
#' @return a 'mareyMap' object with the estimated recombination map.
#' @export
#'
#' @import pbmcapply
#' @import parallel
#' @import utils
#' @import stats
#'
recombinationMap = function(x,
                            chromosome = character(),
                            method = c("loess", "spline"),
                            K = 5,
                            boot = 1000,
                            nCores = 1,
                            smoothing = numeric(),
                            nResampling = 1000,
                            calibrationRange = c(0.2, 0.5),
                            degree = 2,
                            windows = 10^5) {
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
  df = df[which(df$vld == TRUE), ]

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
  x$windows[chromosome] = list(windows)
  x$nBootstrap[chromosome] = boot

  cat("Fit the model.\n")
  if (is.numeric(windows)) {
    boundaries = seq(0, max(df$phys, na.rm = TRUE), by = windows)
    physWindows = data.frame(start = boundaries[-length(boundaries)],
                             end = (boundaries[-1] - 1)) # A sequence of windows coordinates
  } else {
    if (is.data.frame(windows)) {
      colnames(windows) = c("start", "end")
      physWindows = data.frame(start = windows$start,
                               end = windows$end) # A sequence of windows coordinates
      } else {
      warning("'windows' must be an integer of windows size or a dataframe of start/end positions.")
    }
  }

  physWindows$point = round((physWindows$end + physWindows$start)/2, digits = 0)
  bootPrediction = matrix(NA, nrow = nrow(physWindows), ncol = boot)
  mareyPredictions = matrix(NA, nrow = nrow(physWindows), ncol = boot)

  # Fit the Marey function to all points
  if (method == "loess") {
    fitMarey = fitLoess(df, span = smoothingParam, degree = degree)
  }
  if (method == "spline") {
    fitMarey = fitSpline(df, spar = smoothingParam)
  }

  # Estimate a Marey function with Confidence Intervals
  pb = txtProgressBar(min = 1, max = boot, initial = 1)
  for (b in 1:boot) {
    setTxtProgressBar(pb, b)

    resampled = df[sample(1:nrow(df), replace = TRUE), ]

    if (method == "loess") {
      fitMarey = fitLoess(resampled, span = smoothingParam, degree = degree)
      predicted.start = predict(fitMarey, newdata = physWindows$start)
      predicted.end = predict(fitMarey, newdata = physWindows$end)
    }
    if (method == "spline") {
      fitMarey = fitSpline(resampled, spar = smoothingParam)
      predicted.start = predict(fitMarey, physWindows$start)$y
      predicted.end = predict(fitMarey, physWindows$end)$y
    }
    bootPrediction[,b] = ((predicted.end - predicted.start) / (physWindows$end - physWindows$start))
    mareyPredictions[,b] = (predicted.end - predicted.start)
  }
  Sys.sleep(1)
  close(pb)

  cat("Estimate the Marey interpolation function with C.I.\n")
  Y = rowMeans(mareyPredictions, na.rm = TRUE)
  Yupper = apply(mareyPredictions, MARGIN = 1, FUN = upperCI)
  Ylower = apply(mareyPredictions, MARGIN = 1, FUN = lowerCI)
  Y[is.na(Y)] = 0
  Yupper[is.na(Yupper)] = 0
  Ylower[is.na(Ylower)] = 0

  mareyCI = data.frame(
    set = unique(df$set),
    map = chromosome,
    physicalPosition = physWindows$point,
    geneticPositioncM = cumsum(Y),
    upperGeneticPositioncM = cumsum(Y) + (Yupper - Y),
    lowerGeneticPositioncM = cumsum(Y) - (Y - Ylower)
  )

  x$mareyCI = x$mareyCI[x$mareyCI$map != chromosome,]
  x$mareyCI = rbind(x$mareyCI, mareyCI)

  cat("Estimate recombination rates.\n")
  # dX = rowMeans(embed(physWindows$point, 2))
  # upperCI = function(x) {
  #   quantile(x, 0.975, na.rm = TRUE)
  # }
  # lowerCI = function(x) {
  #   quantile(x, 0.025, na.rm = TRUE)
  # }
  dY = rowMeans(bootPrediction, na.rm = TRUE)
  dYupper = apply(bootPrediction, MARGIN = 1, FUN = upperCI)
  dYlower = apply(bootPrediction, MARGIN = 1, FUN = lowerCI)
  dY[dY < 0] = 0
  dYupper[dYupper < 0] = 0
  dYlower[dYlower < 0] = 0

  estimates = data.frame(
    set = unique(df$set),
    map = chromosome,
    start = physWindows$start,
    end = physWindows$end,
    recRate = dY,
    upperRecRate = dYupper,
    lowerRecRate = dYlower
  )

  x$recMap = x$recMap[x$recMap$map != chromosome,]
  x$recMap = rbind(x$recMap, estimates)

  return(x)
}



# Generic functions for Confidence Intervals
upperCI = function(x) {
  quantile(x, 0.975, na.rm = TRUE)
}
lowerCI = function(x) {
  quantile(x, 0.025, na.rm = TRUE)
}
