#' Plot the Lorenz curve
#'
#' @description
#' The Lorenz curve is a graphical representation of the heterogeneity of rates in a dataset, here recombination rates.
#' It plots the ordered distribution of relative recombination rates as a function of their relative physical size.
#' A Lorenz curve always starts at (0,0) and ends at (1,1).
#' The more the curve is close to the diagonal, the less heterogeneous is the distribution of recombination rates.
#'
#'
#' @param x a Marey map object (the Marey map is used).
#' @param return.plot TRUE or FALSE, whether to return a plot or a data frame (default = TRUE returns the plot)
#'
#' @return a data frame with relative physical and genetic position ordered.
#' @export
#'
#'
lorenz = function(x, return.plot = TRUE) {
  df = x$recMap


  df = df[!is.na(df$recRate),]

  df$width = df$end - df$start

  df$relativePhys = df$width/sum(df$width, na.rm = T)
  df$relativeGen = (df$recRate * df$width)/sum((df$recRate * df$width), na.rm = T)

  df = df[order(df$relativeGen),]

  df$relativePhys = cumsum(df$relativePhys)
  df$relativeGen = cumsum(df$relativeGen)

  out = data.frame(set = df$set,
                   map = df$map,
                   relativeGen = df$relativeGen,
                   relativePhys = df$relativePhys)

  diagonal = data.frame(x = seq(0, 1, by = 0.01),
                        y = seq(0, 1, by = 0.01))
  p = ggplot(data = out, aes(x = relativePhys, y = relativeGen)) +
        geom_line() +
        geom_line(data = diagonal, aes(x = x, y = y), color = "Grey") +
        xlim(0, 1) +
        ylim(0, 1) +
        xlab("Proportion of genomic distance") +
        ylab("Proportion of genetic distance") +
        theme_bw()

  if (return.plot) {
    return(p)
  } else {
    return(out)
  }
}
