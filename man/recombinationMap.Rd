% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recombinationMap.R
\name{recombinationMap}
\alias{recombinationMap}
\title{Global function estimate the recombination map}
\usage{
recombinationMap(
  x,
  method = "loess",
  K = 5,
  boot = 1000,
  nCores = 1,
  smoothing = numeric(),
  nResampling = 1000,
  calibrationRange = c(0.2, 0.5),
  degree = 2,
  windows = 10^5,
  setNegativeValues = 0,
  verbose = TRUE
)
}
\arguments{
\item{x}{a 'mareyMap' object.}

\item{method}{an interpolation method, either 'loess' or 'splines' (default = 'loess').}

\item{K}{number of clusters to subset in K-fold cross-validation.}

\item{boot}{number of bootstraps to estimate the confidence interval.}

\item{nCores}{number of cores to parallelize.}

\item{smoothing}{(optional) smoothing parameter, if you do not want to calibrate it.}

\item{nResampling}{number of iterations in the cross-validation procedure.}

\item{calibrationRange}{range of the smoothing parameter space to explore.}

\item{degree}{(optional) the degree parameter of the polynomial used in the 'loess' method.}

\item{windows}{size or coordinates of the windows along the chromosome: either integer in basepairs (default = 10^5) or a data.frame of start/end coordinates.}

\item{setNegativeValues}{the value (e.g. 0 or NA) to use for negative values in recombination rates (default = 0).}

\item{verbose}{Whether to print messages and progress bars at each step (default = TRUE).}
}
\value{
a 'mareyMap' object with the estimated recombination map.
}
\description{
Global function estimate the recombination map
}
