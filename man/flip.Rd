% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flip.R
\name{flip}
\alias{flip}
\title{Flip genetic positions}
\usage{
flip(x)
}
\arguments{
\item{x}{a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame}
}
\value{
the input data frame with flipped genetic positions.
}
\description{
Sometimes the genetic positions and physical positions are not in the same order
(i.e. the map is decreasing). This function allows to flip the genetic positions to get them in the same order as the reference genome
}
