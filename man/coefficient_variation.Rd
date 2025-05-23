% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefficientVariation.R
\name{coefficient_variation}
\alias{coefficient_variation}
\title{Calculate the variance and coefficient of variation of recombination rates}
\usage{
coefficient_variation(x)
}
\arguments{
\item{x}{a `marey_map` object, with a recombination map already estimated}
}
\value{
a list of numeric values of the variance, standard deviation and coefficient of variation.
}
\description{
The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean
CV = standard deviation/mean. It is estimated on recombination rates.
}
