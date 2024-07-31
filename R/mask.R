#' Mask markers in interval
#'
#' @description Mask genetic markers in a Marey map (set vld = FALSE)
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#' @param intervals a data frame of genomic intervals to mask, with physical positions in columns 'start' and 'end'.
#'
#' @return the input data frame with masked genetic positions.
#' @export
#'
mask.interval = function(x, intervals) {
  start = intervals$start
  end = intervals$end
  nintervals = nrow(intervals)

  for (i in 1:nintervals) {
    # Sample markers in the interval
    idx = which(x$mareyMap$phys >= start[i] & x$mareyMap$phys <= end[i])

    if (length(idx) > 0) {
      x$mareyMap[idx, "vld"] = FALSE
    }
  }

  return(x)
}



#' Un-Mask markers in interval
#'
#' @description Mask genetic markers in a Marey map (set vld = TRUE)
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#' @param intervals a data frame of genomic intervals to mask, with physical positions in columns 'start' and 'end'.
#'
#' @return the input data frame with masked genetic positions.
#' @export
#'
unmask.interval = function(x, intervals) {
  start = intervals$start
  end = intervals$end
  nintervals = nrow(intervals)

  for (i in 1:nintervals) {
    # Sample markers in the interval
    idx = which(x$mareyMap$phys >= start[i] & x$mareyMap$phys <= end[i])

    if (length(idx) > 0) {
      x$mareyMap[idx, "vld"] = TRUE
    }
  }

  return(x)
}
