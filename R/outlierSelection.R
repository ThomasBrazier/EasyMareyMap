#' Remove outlying markers by visual zone selection on the Marey map
#'
#' This function allows to graphically remove/keep markers by zone selection.
#'
#' Select a zone by left-click around markers to remove/keep, and close the zone by a right-click.
#'
#' The 'vld' status of the points selected in the graphical windows is set to the opposite:
#' - 'FALSE' if actual status is 'TRUE'
#' - 'TRUE' otherwise
#'
#' @param x a `marey_map` object.
#'
#' @return a `marey_map` object.
#' @export
#'
#' @import gatepoints
#' @import grDevices
#'
outlier_selection = function(x = marey_map()) {
  tmp = x$mareyMap
  idx = point_selection(tmp)
  status = !(x$mareyMap$vld[as.character(row.names(x$mareyMap)) %in% idx])
  x$mareyMap$vld[as.character(row.names(x$mareyMap)) %in% idx] = status
  return(x)
}


point_selection = function(x) {
  df = x[,c("phys", "gen")]
  df = df[!is.na(df$phys),]
  dev.new()
  plot(df, col = "black")
  selectedPoints = gatepoints::fhs(df, mark = TRUE)
  dev.off()
  pointNames = as.character(selectedPoints)
  return(pointNames)
}
