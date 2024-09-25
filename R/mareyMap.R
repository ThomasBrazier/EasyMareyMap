#' Instantiate a `marey_map` object
#'
#' @param x a data frame of a Marey map with columns "set",	"map",	"mkr",	"phys",	"gen",	"vld". If 'chromosome' names is not provided, this data.frame must contain only a single chromosome.
#' @param chromosome the name of the chromosome/map to import from 'x'.
#' @param chromosome_length an optional vector of chromosome length in Mb.
#'
#' @return a `marey_map` object.
#' @export
#'
marey_map = function(x = data.frame(), chromosome = NA, chromosome_length = numeric()) {
  x = as.data.frame(x)
  if (is.na(chromosome)) {
    x = x
  } else {
    if (length(chromosome) == 1) {
      x = x[which(x$map == chromosome),]
    } else {
      if (length(unique(x$map)) != 1) {
        stop(paste0("Only a single chromosome (one map) is admitted.\n"))
      }
    }
  }

  marey = new_marey_map(x, chromosome_length = chromosome_length)
  validate_marey_map(marey)
}

new_marey_map = function(x = data.frame(), chromosome_length = numeric()) {
  stopifnot(is.data.frame(x))
  stopifnot(nrow(x) > 0)
  stopifnot(is.numeric(chromosome_length))
  if (!identical(colnames(x), c("set",	"map",	"mkr",	"phys",	"gen",	"vld"))) {
    stop("Invalid column names in 'x'.\n", paste(colnames(x), sep = ", "))
  }
  if (length(unique(x$map)) != 1) {
    stop(paste0("More than one chromosome in data frame :", paste(unique(x$map), collapse = ', '), "\nOnly a single chromosome (one map) is admitted.\n"))
  }
  x$set = as.factor(x$set)
  x$map = as.factor(x$map)
  x$mkr = as.character(x$mkr)
  x$phys = as.integer(x$phys)
  x$gen = as.integer(x$gen)
  x$vld = as.logical(x$vld)
  x$predict.se = NA

  if (length(chromosome_length) == 0) {
    chrLength = aggregate(x$phys, by = list(x$map), function(x) {max(x, na.rm = TRUE)})
    chromosome_length = chrLength$x
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

  mareyCI = data.frame(set = factor(),
                      map = factor(),
                      physicalPosition = numeric(),
                      geneticPositioncM = numeric(),
                      upperGeneticPositioncM = numeric(),
                      lowerGeneticPositioncM = numeric())

  df = list(mareyMap = x,
            recMap = recMap,
            mareyCI = mareyCI,
            chromosomeName = unique(x$map),
            chromosomeLength = chromosome_length,
            linkageMapLength = linkageMapLength,
            interpolationMethod = character(),
            model = NULL,
            fitted = numeric(),
            residuals = numeric(),
            goodness.r.squared = numeric(),
            goodness.adj.r.squared = numeric(),
            goodness.p = numeric(),
            windows = numeric(),
            smoothingParam = numeric(),
            nBootstrap = integer())

  structure(df,
            class = "marey_map"
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

