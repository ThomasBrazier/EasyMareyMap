#' Global function estimate the recombination map
#'
#' @param x a 'mareyMap' object.
#' @param method an interpolation method, either 'loess' or 'splines' (default = 'loess').
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
                            method = "loess",
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

  nCores = min(nCores,
               (parallel::detectCores(all.tests = FALSE, logical = TRUE)-1))

  cat("Processing...\n")
  df = x$mareyMap
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

  x$interpolationMethod = method
  x$smoothingParam = smoothingParam
  x$windows = list(windows)
  x$nBootstrap = boot

  cat("Fit the model in sliding windows.\n")
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


  # Fit the Marey function to all points
  if (class(fitMarey) == "loess") {
    fitMarey = fitLoess(df, span = smoothingParam, degree = degree)
    x$fitted = fitMarey$fitted
    x$residuals = fitMarey$residuals
    x$model = fitMarey
  }
  if (class(fitMarey) == "smooth.spline") {
    fitMarey = fitSpline(df, spar = smoothingParam)
    x$fitted = fitMarey$y
    x$model = fitMarey
  }



  cat("Estimate a bootstrapped Marey map.\n")
  # Estimate a Marey function with Confidence Intervals
  x$mareyCI = bootstrapMareyMap(x, physWindows, nboot = boot)

  cat("Estimate recombination rates.\n")
  x$recMap = bootstrapRecMap(x, physWindows, nboot = boot)

  return(x)
}



# Generic functions for Confidence Intervals
upperCI = function(x) {
  quantile(x, 0.975, na.rm = TRUE)
}
lowerCI = function(x) {
  quantile(x, 0.025, na.rm = TRUE)
}

#' Bootstrap the Marey map
#'
#' @description A generic function to bootstrap a C.I. for a Marey map on given intervals
#'
#' @param x a 'mareyMap' object.
#' @param intervals a data.frame of genomic windows in which the recombination rate will be estimated, with columns 'start'-'end'.
#' @param nboot number of bootstraps to estimate the confidence interval.
#'
#' @return a data frame of a Marey map with bootstrapped genetic positions.
#' @export
#'
bootstrapMareyMap = function(x, intervals, nboot = boot) {
  stopifnot(class(x) == "mareyMap")

  df = x$mareyMap
  fitMarey = x$model

  mareyPredictions = matrix(NA, nrow = nrow(physWindows), ncol = nboot)

  pb = txtProgressBar(min = 1, max = nboot, initial = 1)
  for (b in 1:boot) {
    setTxtProgressBar(pb, b)

    resampled = df[sample(1:nrow(df), replace = TRUE), ]

    if (class(fitMarey) == "loess") {
      fitMarey = fitLoess(resampled, span = fitMarey$pars$span, degree = fitMarey$pars$degree)
      predicted.start = predict(fitMarey, newdata = physWindows$start)
      predicted.end = predict(fitMarey, newdata = physWindows$end)
    }
    if (class(fitMarey) == "smooth.spline") {
      fitMarey = fitSpline(resampled, spar = fitMarey$spar)
      predicted.start = predict(fitMarey, physWindows$start)$y
      predicted.end = predict(fitMarey, physWindows$end)$y
    }
    mareyPredictions[,b] = (predicted.end - predicted.start)
  }
  Sys.sleep(1)
  close(pb)

  Y = rowMeans(mareyPredictions, na.rm = TRUE)
  Yupper = apply(mareyPredictions, MARGIN = 1, FUN = upperCI)
  Ylower = apply(mareyPredictions, MARGIN = 1, FUN = lowerCI)
  Y[is.na(Y)] = 0
  Yupper[is.na(Yupper)] = 0
  Ylower[is.na(Ylower)] = 0

  mareyCI = data.frame(
    set = unique(df$set),
    map = unique(df$map),
    physicalPosition = physWindows$point,
    geneticPositioncM = cumsum(Y),
    upperGeneticPositioncM = cumsum(Y) + (Yupper - Y),
    lowerGeneticPositioncM = cumsum(Y) - (Y - Ylower)
  )

  return(mareyCI)
}



#' Bootstrap the recombination map
#'
#' @param x a 'mareyMap' object.
#' @param intervals a data.frame of genomic windows in which the recombination rate will be estimated, with columns 'start'-'end'.
#' @param nboot number of bootstraps to estimate the confidence interval.
#'
#' @return a data frame of a recombination map.
#' @export
#'
bootstrapRecMap = function(x, intervals, nboot = boot) {
  stopifnot(class(x) == "mareyMap")

  df = x$mareyMap
  fitMarey = x$model

  bootPrediction = matrix(NA, nrow = nrow(physWindows), ncol = nboot)

  pb = txtProgressBar(min = 1, max = nboot, initial = 1)
  for (b in 1:boot) {
    setTxtProgressBar(pb, b)

    resampled = df[sample(1:nrow(df), replace = TRUE), ]

    if (class(fitMarey) == "loess") {
      fitMarey = fitLoess(resampled, span = fitMarey$pars$span, degree = fitMarey$pars$degree)
      predicted.start = predict(fitMarey, newdata = physWindows$start)
      predicted.end = predict(fitMarey, newdata = physWindows$end)
    }
    if (class(fitMarey) == "smooth.spline") {
      fitMarey = fitSpline(resampled, spar = fitMarey$spar)
      predicted.start = predict(fitMarey, physWindows$start)$y
      predicted.end = predict(fitMarey, physWindows$end)$y
    }
    bootPrediction[,b] = ((predicted.end - predicted.start) / (physWindows$end - physWindows$start))
  }
  Sys.sleep(1)
  close(pb)

  dY = rowMeans(bootPrediction, na.rm = TRUE)
  dYupper = apply(bootPrediction, MARGIN = 1, FUN = upperCI)
  dYlower = apply(bootPrediction, MARGIN = 1, FUN = lowerCI)
  dY[dY < 0] = 0
  dYupper[dYupper < 0] = 0
  dYlower[dYlower < 0] = 0

  estimates = data.frame(
    set = unique(df$set),
    map = unique(df$map),
    start = physWindows$start,
    end = physWindows$end,
    recRate = dY,
    upperRecRate = dYupper,
    lowerRecRate = dYlower
  )

  return(estimates)
}

