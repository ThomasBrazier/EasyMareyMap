#' Remove outlying markers by visual zone selection on the Marey map
#'
#' @param x a mareyMap object.
#' @param chromosome a single chromosome name to display.
#'
#' @return a mareyMap object.
#' @export
#'
#' @import gatepoints
#' @import grDevices
#'
outlierSelection = function(x = mareyMap(), chromosome = character()) {
  tmp = x$mareyMap[which(x$mareyMap$map == chromosome),]
  idx = pointSelection(tmp)
  x$mareyMap$vld[as.character(row.names(x$mareyMap)) %in% idx] = FALSE
  return(x)
}


pointSelection = function(x) {
  df = x[,c("phys", "gen")]
  df = df[!is.na(df$phys),]
  dev.new()
  plot(df, col = "black")
  selectedPoints = gatepoints::fhs(df, mark = TRUE)
  dev.off()
  pointNames = as.character(selectedPoints)
  return(pointNames)
}
