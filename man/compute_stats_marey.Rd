% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparative_marey_map.R
\name{compute_stats_marey}
\alias{compute_stats_marey}
\title{Compute summary statistics on a meta `comparative_marey_map` object}
\usage{
compute_stats_marey(x, statistics = c("mean", "median"), ...)
}
\arguments{
\item{x}{a `comparative_marey_map` object}

\item{statistics}{a vector of statistics to compute}

\item{...}{additional arguments}
}
\value{
a list of summary statistics
}
\description{
Compute summary statistics on a meta `comparative_marey_map` object
}
\details{
The complete list of statistics that can be computed is
`mean`, `median`, `weighted.mean`, `variance`, `gini`, `peripherybias` and `coefficient_variation`
}
\section{Slots}{

\describe{
\item{\code{set}}{a vector of dataset names}

\item{\code{map}}{a vector of chromosome names}

\item{\code{statistics}}{one slot per statistic computed, see details for the complete list of statistics}
}}

