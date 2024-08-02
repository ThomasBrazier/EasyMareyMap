#' Mask markers in interval/position
#'
#' @description Mask genetic markers in a Marey map (set vld = FALSE)
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#' @param position a data frame of genomic intervals/positions to mask, with physical positions in columns 'start' and 'end', or marker positions with columns 'phys' and 'gen'.
#'
#' @return the input data frame with masked genetic positions.
#' @export
#'
mask.marker = function(x, position) {
  x = vld.marker(x, position, vld = FALSE)
  return(x)
}



#' Un-Mask markers in interval/position
#'
#' @description Mask genetic markers in a Marey map (set vld = TRUE)
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#' @param position a data frame of genomic intervals/positions to mask, with physical positions in columns 'start' and 'end', or marker positions with columns 'phys' and 'gen'.
#'
#' @return the input data frame with masked genetic positions.
#' @export
#'
unmask.marker = function(x, position) {
  x = vld.marker(x, position, vld = TRUE)
  return(x)
}





#' Change the state 'vld' of markers in interval/position
#'
#' @description Change 'vld' for genetic markers in a Marey map (set vld = FALSE)
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#' @param position a data frame of genomic intervals/positions to mask, with physical positions in columns 'start' and 'end', or marker positions with columns 'phys' and 'gen'.
#' @param vld the new state for 'vld'.
#'
#' @return the input data frame with new state of 'vld' at markers positions.
#'
vld.marker = function(x, position, vld = FALSE) {
  nposition = nrow(position)

  if (sum(colnames(position) == c("start", "end")) == 2) {
    start = position$start
    end = position$end
    for (i in 1:nposition) {
      # Sample markers in the interval
      idx = which(x$mareyMap$phys >= start[i] & x$mareyMap$phys <= end[i])
      if (length(idx) > 0) {
        x$mareyMap[idx, "vld"] = FALSE
      }
    }
  } else {
    if (sum(colnames(position) == c("phys", "gen")) == 2) {
      phys = position$phys
      gen = position$gen
      for (i in 1:nposition) {
        # Sample markers in the interval
        idx = which(x$mareyMap$phys == phys[i] & x$mareyMap$gen == gen[i])
        if (length(idx) > 0) {
          x$mareyMap[idx, "vld"] = FALSE
        }
      }
    }
  }
  return(x)
}
