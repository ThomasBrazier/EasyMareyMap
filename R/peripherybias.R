#' Calculate the periphery-bias ratio as in Brazier and Glémin (2022)
#'
#' @param x a `marey_map` object, with a recombination map already estimated.
#' @param proportion_periphery the proportion of the periphery to keep (default = 0.1)
#' @param chromosome_arm which chromosome arm (half of the chromosome) to keep (either 'random', 'left', 'right' or 'both')
#'
#' @return a numeric value of the periphery-bias ratio.
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#' 
peripherybias = function(x, proportion_periphery = 0.1,  chromosome_arm = "random") {
  df = x$recMap
  chrlength = x$chromosomeLength

  # Cut in one hundred segments with their recombination rate
  intervals = GRanges(seqnames = x$chromosomeName,
                      ranges = IRanges(start = seq(1, chrlength, round(chrlength/100, digits = 0)), width = round(chrlength/100, digits = 0)),
                      strand = "*")

  hundred_windows = windows_physical_map(x, windows = intervals)

  # Keep only one side (one chromosome arm) or both
  # Mask with NAs

  stopifnot(chromosome_arm %in% c("random", "left", "right", "both"))

  if (chromosome_arm == "random") {
    chromosome_arm = sample(c("left", "right"), 1)
  }

  if (chromosome_arm == "left") {
    mask = c(51:100)
  } else {
    if (chromosome_arm == "right") {
      mask = c(1:50)
    }
  }

  if (chromosome_arm != "both") {
    hundred_windows$recRate[mask] = NA
    hundred_windows$upperRecRate[mask] = NA
    hundred_windows$lowerRecRate[mask] = NA
  }

  # The periphery-bias ratio is the mean recombination rate
  # of the proportion sampled at the periphery divided by the mean of the whole chromosome arm
  periphery_region = c(seq(1, round(length(intervals) * proportion_periphery, digits = 0), by = 1),
                       seq(round(100 - length(intervals) * proportion_periphery, digits = 0), 100, by = 1))

  peripherybias_ratio = (mean(hundred_windows$recRate[periphery_region], na.rm = TRUE)/mean(hundred_windows$recRate, na.rm = TRUE))

  return(peripherybias_ratio)
}
