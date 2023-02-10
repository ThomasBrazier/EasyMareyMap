#' Remove outlying markers by visual zone selection on the Marey map
#'
#' @param x A mareyMap object
#' @param chromosome A single chromosome name to display
#'
#' @return A mareyMap object
#' @export
#'
#' @examples
outlierSelection = function(x = mareyMap(), chromosome = character()) {
  tmp = x$mareyMap[which(x$mareyMap$map == chromosome),]
  idx = pointSelection(tmp)
  x$mareyMap$vld[as.character(row.names(x$mareyMap)) %in% idx] = FALSE
  return(x)
}


pointSelection = function(x) {
  df = x[,c("phys", "gen")]
  df = df[!is.na(df$phys),]
  X11()
  plot(df, col = "black")
  selectedPoints = fhs(df, mark = TRUE)
  dev.off()
  pointNames = as.character(selectedPoints)
  return(pointNames)
}
