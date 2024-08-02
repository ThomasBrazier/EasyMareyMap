#' Aggregate recombination rates in sliding windows along the genome
#'
#' @param x a Marey map object, with a recombination map already estimated.
#' @param windows a GRanges object with genomic intervals
#' @param method the method used to aggregate (either "mean", "weightedmean" or "median")
#'
#' @return a new data frame with recombination rates in sliding windows.
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
windowsPhysicalMap = function(x, windows = GRanges(), method = "mean") {
  df = x$recMap

  recMap = makeGRangesFromDataFrame(df,
                                    seqnames.field = "map",
                                    keep.extra.columns = T)

  hits = findOverlaps(windows, recMap)

  windows$set = unique(df$set)

  # Aggregate recombination rates in intervals
  windows$recRate = NA
  windows$upperRecRate = NA
  windows$lowerRecRate = NA

  if (method == "mean") {
    windows$recRate = unlist(lapply(1:length(windows), function(x) {mean(recMap$recRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
    windows$upperRecRate = unlist(lapply(1:length(windows), function(x) {mean(recMap$upperRecRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
    windows$lowerRecRate = unlist(lapply(1:length(windows), function(x) {mean(recMap$lowerRecRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
  }
  if (method == "weightedmean") {
    windows$recRate = unlist(lapply(1:length(windows), function(x) {
      overlaps = recMap[subjectHits(hits)[which(queryHits(hits) == x)]]
      overlaps = unlist(as(lapply(1:length(overlaps),
                                  function(y) {restrict(overlaps[y],
                                                        start = start(intervals)[x],
                                                        end = end(intervals)[x])}),
                           "GRangesList"))
      return(weighted.mean(overlaps$recRate, width(overlaps), na.rm = TRUE))}))

    windows$upperRecRate = unlist(lapply(1:length(windows), function(x) {
      overlaps = recMap[subjectHits(hits)[which(queryHits(hits) == x)]]
      overlaps = unlist(as(lapply(1:length(overlaps),
                                  function(y) {restrict(overlaps[y],
                                                        start = start(intervals)[x],
                                                        end = end(intervals)[x])}),
                           "GRangesList"))
      return(weighted.mean(overlaps$upperRecRate, width(overlaps), na.rm = TRUE))}))

    windows$lowerRecRate = unlist(lapply(1:length(windows), function(x) {
      overlaps = recMap[subjectHits(hits)[which(queryHits(hits) == x)]]
      overlaps = unlist(as(lapply(1:length(overlaps),
                                  function(y) {restrict(overlaps[y],
                                                        start = start(intervals)[x],
                                                        end = end(intervals)[x])}),
                           "GRangesList"))
      return(weighted.mean(overlaps$lowerRecRate, width(overlaps), na.rm = TRUE))}))
  }
  if (method == "median") {
    windows$recRate = unlist(lapply(1:length(windows), function(x) {median(recMap$recRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
    windows$upperRecRate = unlist(lapply(1:length(windows), function(x) {median(recMap$upperRecRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
    windows$lowerRecRate = unlist(lapply(1:length(windows), function(x) {median(recMap$lowerRecRate[subjectHits(hits)[which(queryHits(hits) == x)]], na.rm = TRUE)}))
  }

  windows = as.data.frame(windows)

  df2 = data.frame(set = windows$set,
                   map = windows$seqnames,
                   start = windows$start,
                   end = windows$end,
                   recRate = windows$recRate,
                   upperRecRate = windows$upperRecRate,
                   lowerRecRate = windows$lowerRecRate)
  return(df2)
}
