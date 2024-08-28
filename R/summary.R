#' Summary of a Marey map object
#'
#' @param object a Marey map object.
#' @param ... arguments passed to the generic summary function.
#'
#' @import ggplot2
#'
#' @return a summary of the Marey map object.
#'
#' @method summary mareyMap
#' @export
summary.mareyMap = function(object, ...) {
  dataset = levels(object$mareyMap$set)
  nameChromosome = object$chromosomeName
  lengthChromosome = object$chromosomeLength/10^6
  lengthLinkageMap = object$linkageMapLength
  nbMarkers = sum(object$mareyMap$vld == TRUE)
  densityMarkers = nbMarkers/sum(lengthChromosome)
  chromEst = names(object$smoothingParam)
  if (nrow(object$recMap) > 0) {
    recRate = aggregate(recRate ~ map, object$recMap, function(x) {summary(x, digits = 3)})
  } else {
    recRate = list()
  }
  smoothingParam = object$smoothingParam
  interpolation = object$interpolationMethod
  nbWindows = nrow(object$recMap)
  nBootstrap = object$nBootstrap
  cat("============== Summary of the Marey map ==============\n",
      "Dataset: ", dataset, "\n",
      "Chromosome name: ", nameChromosome, "\n",
      "Chromosome length (Mb): ", lengthChromosome, "\n",
      "Linkage map length (cM): ", lengthLinkageMap, "\n",
      "Number of markers: ", nbMarkers, "\n",
      "Markers density (marker/Mb): ", densityMarkers, "\n",
      "============== Recombination map ==============\n",
      "Recombination rates (cM/Mb):\n",
      "Min.   1st Qu.    Median      Mean   3rd Qu.      Max.\n")
  if (length(recRate) > 0) {for (i in 1:nrow(recRate)) {cat(unlist(recRate[i,]) * 1e6, "\n")}}
  cat("============== Interpolation ==============\n")
  cat("Interpolation method: ", interpolation, "\n",
      "Smoothing parameter: ", smoothingParam, "\n",
      "Number of windows: ", nbWindows, "\n",
      "Number of bootstraps: ", nBootstrap, "\n")
  cat("============== Goodness of fit (fitted ~ y) ==============\n",
      "fitted ~ y R squared: ", object$goodness.r.squared, "\n",
      "fitted ~ y adjusted R squared: ", object$goodness.adj.r.squared, "\n",
      "fitted ~ y p-value ", object$goodness.p, "\n")
}
