#' Global function estimate the recombination map
#'
#' @param x a `marey_map` object.
#' @param method an interpolation method, either 'loess' or 'splines' (default = 'loess').
#' @param K number of clusters to subset in K-fold cross-validation.
#' @param boot number of bootstraps to estimate the confidence interval.
#' @param n_cores number of cores to parallelize.
#' @param smoothing (optional) smoothing parameter, if you do not want to calibrate it.
#' @param n_resampling number of iterations in the cross-validation procedure.
#' @param calibration_range range of the smoothing parameter space to explore.
#' @param degree (optional) the degree parameter of the polynomial used in the 'loess' method.
#' @param windows size or coordinates of the windows along the chromosome: either integer in basepairs (default = 10^5) or a data.frame of start/end coordinates.
#' @param set_negative_values the value (e.g. 0 or NA) to use for negative values in recombination rates (default = 0).
#' @param verbose Whether to print messages and progress bars at each step (default = TRUE).
#'
#' @return a `marey_map` object with the estimated recombination map.
#' @export
#'
#' @import pbmcapply
#' @import parallel
#'
recombination_map = function(x,
                            method = "loess",
                            K = 5,
                            boot = 1000,
                            n_cores = 1,
                            smoothing = numeric(),
                            n_resampling = 1000,
                            calibration_range = c(0.2, 0.5),
                            degree = 2,
                            windows = 10^5,
                            set_negative_values = 0,
                            verbose = TRUE) {
  if (!(method %in% c("loess", "spline"))) {
    stop("Interpolation method must be specified: loess or spline.")
  }

  calibrationFrom = calibration_range[1]
  calibrationTo = calibration_range[2]

  nCores = min(n_cores,
               (parallel::detectCores(all.tests = FALSE, logical = TRUE)-1))

  if(verbose) {cat("Processing...\n")}
  df = get_marey_map(x)

  if (length(smoothing) > 0) {
    if (is.numeric(smoothing)) {
      smoothingParam = smoothing
    }
  } else {
    if(verbose) {cat("Calibrate the smoothing parameter.\n")}
    smoothingParam = calibrate_smoothing(df,
                                        method = method,
                                        n_resampling = n_resampling,
                                        K = K,
                                        from = calibrationFrom,
                                        to = calibrationTo,
                                        n_cores = nCores,
                                        degree = degree,
                                        verbose
    )
  }

  x$interpolationMethod = method
  x$smoothingParam = smoothingParam
  x$windows = list(windows)
  x$nBootstrap = boot

  if(verbose) {cat("Fit the model in sliding windows.\n")}
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
  if (x$interpolationMethod == "loess") {
    fitMarey = fit_loess(df, span = smoothingParam, degree = degree)
    x$fitted = fitMarey$fitted
    x$residuals = fitMarey$residuals
    x$model = fitMarey
  }
  if (x$interpolationMethod == "spline") {
    fitMarey = fit_spline(df, spar = smoothingParam)
    x$fitted = fitMarey$y
    x$residuals = residuals(fitMarey)
    x$model = fitMarey
  }
  # Goodness of fit criterion
  # OLS regression fitted genetic positions ~ true genetic positions
  goodnessFit = lm(x$fitted ~ x$mareyMap$phys)
  summary(goodnessFit)

  x$goodness.r.squared = as.numeric(summary(goodnessFit)["r.squared"])
  x$goodness.adj.r.squared = as.numeric(summary(goodnessFit)["adj.r.squared"])
  x$goodness.p = as.numeric(unlist(anova(goodnessFit)["Pr(>F)"])[1])


  if(verbose) {cat("Estimate a bootstrapped Marey map.\n")}
  # Estimate a Marey function with Confidence Intervals
  x$mareyCI = bootstrap_marey_map(x, physWindows, nboot = boot, verbose)

  if(verbose) {cat("Estimate recombination rates.\n")}
  x$recMap = bootstrap_rec_map(x, physWindows, nboot = boot, set_negative_values, verbose)


  return(x)
}



# Generic functions for Confidence Intervals
#' @importFrom stats quantile
upperCI = function(x) {
  quantile(x, 0.975, na.rm = TRUE)
}

#' @importFrom stats quantile
lowerCI = function(x) {
  quantile(x, 0.025, na.rm = TRUE)
}

#' Bootstrap the Marey map
#'
#' @description A generic function to bootstrap a C.I. for a Marey map on given intervals
#'
#' @param x a `marey_map` object.
#' @param intervals a data.frame of genomic windows in which the recombination rate will be estimated, with columns 'start'-'end'.
#' @param nboot number of bootstraps to estimate the confidence interval.
#' @param verbose Whether to print messages and progress bars at each step (default = TRUE).
#'
#' @return a data frame of a Marey map with bootstrapped genetic positions.
#' @export
#'
#' @importFrom methods is
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats predict
#' @importFrom stats residuals
#' @importFrom stats anova
#'
bootstrap_marey_map = function(x, intervals, nboot = 1000, verbose = TRUE) {
  stopifnot(is(x, "marey_map"))

  df = x$mareyMap
  fitMarey = x$model

  mareyPredictions = matrix(NA, nrow = nrow(intervals), ncol = nboot)

  if(verbose) {pb = txtProgressBar(min = 1, max = nboot, initial = 1)}
  for (b in 1:nboot) {
    if(verbose) {setTxtProgressBar(pb, b)}

    resampled = df[sample(1:nrow(df), replace = TRUE), ]

    if (is(fitMarey, "loess")) {
      fitMarey = fit_loess(resampled, span = fitMarey$pars$span, degree = fitMarey$pars$degree)
      predicted.start = predict(fitMarey, newdata = intervals$start)
      predicted.end = predict(fitMarey, newdata = intervals$end)
    }
    if (is(fitMarey, "ss")) {
      fitMarey = fit_spline(resampled, spar = fitMarey$spar)
      predicted.start = predict(fitMarey, intervals$start)$y
      predicted.end = predict(fitMarey, intervals$end)$y
    }
    mareyPredictions[,b] = (predicted.end - predicted.start)
  }
  Sys.sleep(1)
  if(verbose) {close(pb)}

  Y = rowMeans(mareyPredictions, na.rm = TRUE)
  Yupper = apply(mareyPredictions, MARGIN = 1, FUN = upperCI)
  Ylower = apply(mareyPredictions, MARGIN = 1, FUN = lowerCI)
  Y[is.na(Y)] = 0
  Yupper[is.na(Yupper)] = 0
  Ylower[is.na(Ylower)] = 0

  mareyCI = data.frame(
    set = unique(df$set),
    map = unique(df$map),
    physicalPosition = intervals$point,
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
#' @param set_negative_values the value (e.g. 0 or NA) to use for negative values in recombination rates (default = 0).
#' @param verbose Whether to print messages and progress bars at each step (default = TRUE).
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats predict
#' 
#' @return a data frame of a recombination map.
#' @export
#'
bootstrap_rec_map = function(x, intervals, nboot = 1000, set_negative_values = 0, verbose = TRUE) {
  stopifnot(class(x) == "marey_map")

  df = x$mareyMap
  fitMarey = x$model

  bootPrediction = matrix(NA, nrow = nrow(intervals), ncol = nboot)

  if(verbose) {pb = txtProgressBar(min = 1, max = nboot, initial = 1)}
  for (b in 1:nboot) {
    if(verbose) {setTxtProgressBar(pb, b)}

    resampled = df[sample(1:nrow(df), replace = TRUE), ]

    if (is(fitMarey, "loess")) {
      fitMarey = fit_loess(resampled, span = fitMarey$pars$span, degree = fitMarey$pars$degree)
      predicted.start = predict(fitMarey, newdata = intervals$start)
      predicted.end = predict(fitMarey, newdata = intervals$end)
    }
    if (is(fitMarey, "ss")) {
      fitMarey = fit_spline(resampled, spar = fitMarey$spar)
      predicted.start = predict(fitMarey, intervals$start)$y
      predicted.end = predict(fitMarey, intervals$end)$y
    }
    bootPrediction[,b] = ((predicted.end - predicted.start) / (intervals$end - intervals$start))
  }
  Sys.sleep(1)
  if(verbose) {close(pb)}

  dY = rowMeans(bootPrediction, na.rm = TRUE)
  dYupper = apply(bootPrediction, MARGIN = 1, FUN = upperCI)
  dYlower = apply(bootPrediction, MARGIN = 1, FUN = lowerCI)
  dY[dY < 0] = set_negative_values
  dYupper[dYupper < 0] = set_negative_values
  dYlower[dYlower < 0] = set_negative_values

  estimates = data.frame(
    set = unique(df$set),
    map = unique(df$map),
    start = intervals$start,
    end = intervals$end,
    recRate = dY,
    upperRecRate = dYupper,
    lowerRecRate = dYlower
  )

  return(estimates)
}

