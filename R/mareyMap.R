#' A Marey map object
#'
#' @param x a data frame of a Marey map with columns "set",	"map",	"mkr",	"phys",	"gen",	"vld".
#' @param chromosomeLength an optional vector of chromosome length in Mb.
#'
#' @return a 'mareyMap' object.
#' @export
#'
mareyMap = function(x = data.frame(), chromosomeLength = numeric()) {
  x = as.data.frame(x)
  marey = new_marey_map(x, chromosomeLength = chromosomeLength)
  validate_marey_map(marey)
}

new_marey_map = function(x = data.frame(), chromosomeLength = numeric()) {
  stopifnot(is.data.frame(x))
  stopifnot(is.numeric(chromosomeLength))
  if (!identical(colnames(x), c("set",	"map",	"mkr",	"phys",	"gen",	"vld"))) {
    stop("Invalid column names in 'x'.\n", colnames(x))
  }
  x$set = as.factor(x$set)
  x$map = as.factor(x$map)
  x$mkr = as.character(x$mkr)
  x$phys = as.integer(x$phys)
  x$gen = as.integer(x$gen)
  x$vld = as.logical(x$vld)

  if (length(chromosomeLength) == 0) {
    chrLength = aggregate(x$phys, by = list(x$map), function(x) {max(x, na.rm = TRUE)})
    chromosomeLength = chrLength$x
  }

  linkageLength = aggregate(x$gen, by = list(x$map), function(x) {max(x, na.rm = TRUE)})
  linkageMapLength = linkageLength$x

  recMap = data.frame(set = factor(),
                      map = factor(),
                      start = integer(),
                      end = integer(),
                      recRate = numeric(),
                      lowerRecRate = numeric(),
                      upperRecRate = numeric())

  df = list(mareyMap = x,
            recMap = recMap,
            chromosomeName = unique(x$map),
            chromosomeLength = chromosomeLength,
            linkageMapLength = linkageMapLength,
            interpolationMethod = character(),
            windows = numeric(),
            smoothingParam = numeric(),
            nBootstrap = integer())

  structure(df,
            class = "mareyMap"
  )
}

validate_marey_map = function(x) {
  values = unclass(x)
  if (length(x$chromosomeName) != length(x$chromosomeLength)) {
    stop("Chromosome name and chromosome length are not of same size.\n",
         paste(x$chromosomeName, collapse = ", "), "\n",
         paste(x$chromosomeLength, collapse = ", "))
  }
  x
}

