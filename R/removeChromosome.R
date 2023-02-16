#' Remove a complete chromosome of the dataset (markers set to vld = FALSE)
#'
#' @param x A mareyMap object
#' @param chromosome A single chromosome name
#'
#' @return A mareyMap object
#' @export
#'
#' @examples removeChromosome(x, "1")
removeChromosome = function(x, chromosome = character()) {
  x$mareyMap$vld[which(x$mareyMap$map == chromosome)] = FALSE
  return(x)
}
