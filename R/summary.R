#' Summary of a Marey map object
#'
#' @param object a Marey map object.
#' @param ... arguments passed to the generic summary function.
#'
#' @return a summary of the Marey map object.
#' @export
#'
#' @import ggplot2
#'
summary.mareyMap = function(object, ...) {
  dataset = levels(object$mareyMap$set)
  nbChromosomes = length(object$chromosomeName)
  nameChromosomes = object$chromosomeName
  lengthChromosomes = object$chromosomeLength/10^6
  lengthLinkageMap = object$linkageMapLength
  nbMarkers = sum(object$mareyMap$vld == TRUE)
  densityMarkers = nbMarkers/sum(lengthChromosomes)
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
      "Number of chromosomes: ", nbChromosomes, "\n",
      "Chromosome names: ", nameChromosomes, "\n",
      "Chromosome length (Mb): ", lengthChromosomes, "\n",
      "Linkage map length (cM): ", lengthLinkageMap, "\n",
      "Number of markers: ", nbMarkers, "\n",
      "Markers density (marker/Mb): ", densityMarkers, "\n",
      "============== Recombination map ==============\n",
      "Chromosomes: ", chromEst, "\n",
      "Recombination rates:\n",
      "Min.   1st Qu.    Median      Mean   3rd Qu.      Max.\n")
  if (length(recRate) > 0) {for (i in 1:nrow(recRate)) {cat(unlist(recRate[i,]), "\n")}}
  cat("Interpolation method: ", interpolation, "\n",
      "Smoothing parameter: ", smoothingParam, "\n",
      "Number of windows: ", nbWindows, "\n",
      "Number of bootstraps: ", nBootstrap, "\n")
}
