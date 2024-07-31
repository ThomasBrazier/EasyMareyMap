#' Flip genetic positions
#'
#' @description Sometimes the genetic positions and physical positions are not in the same order
#' (i.e. the map is decreasing). This function allows to flip the genetic positions to get them in the same order as the reference genome
#'
#' @param x a data frame of a given interval/selection/chromosome with same columns as a 'mareyMap' data frame
#'
#' @return the input data frame with flipped genetic positions.
#' @export
#'
flip = function(x) {
  # A vector of genetic distances values to work on
  vec = x$gen
  # Flip the selection selected segment
  new_gendist = abs(vec - max(vec, na.rm = TRUE)) + min(vec, na.rm = TRUE)
  # Reaffect values to the map
  x$gen = new_gendist
  return(x)
}


#' Flip genetic positions in interval
#'
#' @description This function allows to flip the genetic positions in a given genomic interval
#' to get them in the same order as the reference genome. The intervals are given in a data frame with columns 'start' and 'end'.
#' Each row is an interval with start and end positions in bp.
#' You can flip the whole chromosome by giving a single interval with an end equal or larger than chromosome size.
#'
#' @usage flip.interval(marey, intervals = data.frame(start = 1000, end = 90000))
#'
#' @param x a 'mareyMap' object.
#' @param intervals a data frame of genomic intervals to flip, with physical positions in columns 'start' and 'end'.
#'
#' @return the input 'mareyMap' object with flipped genetic positions.
#' @export
#'
flip.interval = function(x, intervals) {
  start = intervals$start
  end = intervals$end
  nintervals = nrow(intervals)

  for (i in 1:nintervals) {
    # Sample markers in the interval
    idx = which(x$mareyMap$phys >= start[i] & x$mareyMap$phys <= end[i])

    if (length(idx) < 2) {
      stop("No genetic positions to flip")
    }

    # Flip them and replace positions in the MareyMap object
    x$mareyMap[idx,] = flip(x$mareyMap[idx,])
  }

  return(x)
}


#' Flip genetic positions with visual selection of markers
#'
#' @description This function allows to flip the genetic positions by visually selecting the segment to flip.
#'
#' @usage flip.selection(marey)
#'
#' @param x a 'mareyMap' object.
#'
#' @return the input 'mareyMap' object with flipped genetic positions.
#' @export
#'
flip.selection = function(x) {
  marey = x$mareyMap
  idx = pointSelection(marey)
  if (length(idx) < 2) {
    stop("No genetic positions to flip")
  }
  # A vector of genetic distances values to work on
  vec = marey$gen[as.numeric(idx)]
  # Flip the selection selected segment
  new_gendist = abs(vec - max(vec, na.rm = TRUE)) + min(vec, na.rm = TRUE)
  # Reaffect values to the map
  marey$gen[as.numeric(idx)] = new_gendist

  x$mareyMap = marey
  return(x)
}
