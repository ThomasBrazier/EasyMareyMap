summary.mareyMap = function(x) {
  dataset = levels(x$mareyMap$set)
  nbChromosomes = length(x$chromosomeName)
  nameChromosomes = x$chromosomeName
  lengthChromosomes = x$chromosomeLength/10^6
  lengthLinkageMap = x$linkageMapLength
  nbMarkers = sum(x$mareyMap$vld == TRUE)
  densityMarkers = nbMarkers/sum(lengthChromosomes)
  chromEst = names(x$smoothingParam)
  if (nrow(x$recMap) > 0) {
    recRate = aggregate(recRate ~ map, x$recMap, function(x) {summary(x, digits = 3)})
  } else {
    recRate = list()
  }
  smoothingParam = x$smoothingParam
  interpolation = x$interpolationMethod
  windowSize = x$windowSize
  nbWindows = nrow(x$recMap)
  nBootstrap = x$nBootstrap
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
      "Windows size: ", windowSize, "\n",
      "Number of windows: ", nbWindows, "\n",
      "Number of bootstraps: ", nBootstrap, "\n")
}
