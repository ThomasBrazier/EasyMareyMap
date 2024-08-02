#' Calculate the r-bar statistic from Veller et al. (2019)
#' r-bar estimate the intra-chromosomal shuffling
#'
#' @param x a Marey map object.
#' @param genomeSize The total genome size must be provided to estiate the intra-chromosomal genetic shuffling.
#' @param mappingFunction The mapping function used to transform genetic distances (either "none", "haldane" or "kosambi")
#'
#' @return a numeric value of r-bar.
#' @export
#'
#'
veller = function(x, genomeSize = NA, mappingFunction = "none") {
  d = getMareyMap(x)
  # Chromosome length
  chrlength = x$chromosomeLength

  # Compute inter-loci distances for each locus pair (i,j) in a matrix
  dij = matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      dij[i,j] = d$gen[j] - d$gen[i]
    }
  }
  # 0 distances in the diagonal passed to NA
  dij[which(dij == 0)] = NA

  # Convert distances back to recombinant fractions
  if (mappingFunction == "none") {
    # reversed Morgan, converting Morgans to recombinant fraction
    cat("Converting Morgans to recombinant fraction\n")
    rij = dij/100
    # cat(mean(rij, na.rm = TRUE))
  } else {
    if (mappingFunction == "haldane") {
      # Reversed Haldane, converting cM to recombinant fraction
      cat("Reversed Haldane, converting cM to recombinant fraction\n")
      rij = 0.5*(1-exp(-2*dij/100))
      # cat(mean(rij, na.rm = TRUE))
    } else {
      if (mappingFunction == "kosambi") {
        # Reversed Kosambi, converting cM to recombinant fraction
        cat("Reversed Kosambi, converting cM to recombinant fraction\n")
        rij = 0.5*tanh(2*dij/100)
        # cat(mean(rij, na.rm = TRUE))
      }
    }
  }

  # Sum of rij for i<j, i.e. the upper triangular matrix rij
  rij[lower.tri(rij)] = NA

    # Intra-chromosomal component is the averaged rate of shuffling rij for each locus pair (i,j)
  # lambda = nrow(d)*(nrow(d)-1)/2 # loci x (loci-1)/2
  lambda = sum(!is.na(rij)) # is the total number of pairs of loci, i.e. number of values not NA
  r_intra = (sum(rij, na.rm = TRUE)/(lambda*(lambda-1)*(1/2)))*(chrlength/genomeSize)^2
  return(r_intra)
}

