% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/windowsPhysicalMap.R
\name{windowsPhysicalMap}
\alias{windowsPhysicalMap}
\title{Aggregate recombination rates in sliding windows along the genome}
\usage{
windowsPhysicalMap(x, windows = GRanges(), method = "mean")
}
\arguments{
\item{x}{a Marey map object, with a recombination map already estimated.}

\item{windows}{a GRanges object with genomic intervals}

\item{method}{the method used to aggregate (either "mean", "weightedmean" or "median")}
}
\value{
a new data frame with recombination rates in sliding windows.
}
\description{
Aggregate recombination rates in sliding windows along the genome
}
