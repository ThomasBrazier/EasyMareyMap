#' Get the Marey map of a mareyMap object
#'
#' @description
#' Get the Marey map slot and remove markers where `vld` == FALSE. Also remove rows with NAs.
#'
#' @param x a Marey map object.
#'
#' @return a Marey map without markers `vld` == FALSE.
getMareyMap = function(x) {
  df = x$mareyMap
  df = df[which(df$vld == TRUE),]
  df = df[!is.na(df$gen), ]
  df = df[!is.na(df$phys), ]
  return(df)
}
