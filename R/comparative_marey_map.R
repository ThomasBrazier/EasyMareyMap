#' A comparative Marey map object
#'
#' @param x a data frame of a Marey map with columns "set",	"map",	"mkr",	"phys",	"gen",	"vld". If 'chromosome' names is not provided, this data.frame must contain only a single chromosome.
#' @param chromosome the name of the chromosome/map to import from 'x'.
#' @param chromosomeLength an optional vector of chromosome length in Mb.
#'
#' @return a 'mareyMap' object.
#' @export
#'
comparative_mareyMap = function(x = data.frame(), chromosome = NA, chromosomeLength = numeric()) {
  
}