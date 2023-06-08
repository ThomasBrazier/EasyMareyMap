#' Remove a complete chromosome of the dataset (markers set to vld = FALSE)
#'
#' @param x a mareyMap object.
#' @param chromosome a single chromosome name.
#'
#' @return a mareyMap object.
#' @export
#'
removeChromosome = function(x, chromosome = character()) {
  x$mareyMap$vld[which(x$mareyMap$map == chromosome)] = FALSE
  return(x)
}
