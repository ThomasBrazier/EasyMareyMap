#' Plot the Lorenz curve
#'
#' @description
#' The Lorenz curve is a graphical representation of the heterogeneity of rates in a dataset, here recombination rates.
#' It plots the ordered distribution of relative recombination rates as a function of their relative physical size.
#' A Lorenz curve always starts at (0,0) and ends at (1,1).
#' The more the curve is close to the diagonal, the less heterogeneous is the distribution of recombination rates.
#'
#'
#' @param x a Marey map object (the Marey map is used) or a `comparative_marey_map` object
#' @param return.plot TRUE or FALSE, whether to return a plot or a data frame (default = TRUE returns the plot)
#'
#' @return a data frame with relative physical and genetic position ordered.
#' @export
#' @importFrom ggplot2 .data
#'
lorenz = function(x, return.plot = TRUE) {
  if (class(x) == 'mareyMap') {
    df = x$recMap
  } else {
    if (class(x) == 'comparative_marey_map') {
      df = comparative_recmap_to_dataframe(x)
    }
  }
  
  df = df[!is.na(df$recRate),]
  
  df$width = df$end - df$start
  
  # Compute total length (sum of distances)
  total_dist = function(x) {sum(df$width[which(df$set == df$set[x] & df$map == df$map[x])], na.rm = TRUE)}
  total_dist_est = unlist(lapply(c(1:nrow(df)), FUN = total_dist))
  df$relativePhys = df$width/total_dist_est
  
  total_gen = function(x) {
    idx = which(df$set == df$set[x] & df$map == df$map[x])
    sum((df$recRate[idx] * df$width[idx]), na.rm = TRUE)
    }
  total_gen_est = unlist(lapply(c(1:nrow(df)), FUN = total_gen))
  df$relativeGen = (df$recRate * df$width)/total_gen_est
  
  df = df[order(df$set, df$map, df$relativeGen),]
  
  for (i in unique(df$set)) {
    for (j in unique(df$map)) {
      idx = which(df$set == i & df$map == j)
      if (length(idx) > 0) {
        df$relativePhys[idx] = cumsum(df$relativePhys[idx])
        df$relativeGen[idx] = cumsum(df$relativeGen[idx])
      }
    }
  }
  
  # df$relativePhys = cumsum(df$relativePhys)
  # df$relativeGen = cumsum(df$relativeGen)
  
  out = data.frame(set = df$set,
                   map = df$map,
                   relativeGen = df$relativeGen,
                   relativePhys = df$relativePhys)
  
  diagonal = data.frame(x = seq(0, 1, by = 0.01),
                        y = seq(0, 1, by = 0.01))
  
  p = ggplot(data = out, aes(x = relativePhys, y = relativeGen, fill = as.factor(map), colour = as.factor(set))) +
    geom_line(alpha = 0.3) +
    # facet_wrap(~ as.factor(set)) +
    geom_line(data = diagonal, aes(x = x, y = y, fill = NA, colour = NA), color = "black") +
    xlim(0, 1) +
    ylim(0, 1) +
    xlab("Proportion of genomic distance") +
    ylab("Proportion of genetic distance") +
    labs(colour = "Set") +
    theme_bw()
  
  if (return.plot) {
    return(p)
  } else {
    return(out)
  }
}
 