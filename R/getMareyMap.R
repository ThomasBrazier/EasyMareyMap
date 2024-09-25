#' Get the Marey map of a `marey_map` object
#'
#' @description
#' Get the Marey map slot and remove markers where `vld` == FALSE. Also remove rows with NAs.
#'
#' @param x a `marey_map` object.
#'
#' @return a Marey map without markers `vld` == FALSE.
get_marey_map = function(x) {
  df = x$mareyMap
  df = df[which(df$vld == TRUE),]
  df = df[!is.na(df$gen), ]
  df = df[!is.na(df$phys), ]
  return(df)
}
