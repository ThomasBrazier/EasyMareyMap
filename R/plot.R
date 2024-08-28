#' Plot the Marey maps
#'
#' @param x a Marey map object.
#' @param return.plot TRUE or FALSE, whether to return a plot or a data frame (default = TRUE returns the plot)
#' @param ... arguments passed to the generic summary function.
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#'
#' @return a plot of the Marey maps, or a data frame.
#'
#' @method plot mareyMap
#' @export
plot.mareyMap = function(x, return.plot = TRUE, ...) {
  df = x$mareyMap
  # Add the Marey interpolated function
  marey = x$mareyCI
  if (nrow(marey) > 0) {
    marey$vld = TRUE

    p = ggplot2::ggplot(data = df, aes(x = .data$phys/10^6, y = .data$gen, colour = .data$vld)) +
      geom_point(alpha = 0.4) +
      scale_colour_manual(values=c("TRUE" = "dodgerblue4", "FALSE" = "firebrick4")) +
      facet_grid(~as.factor(set)) +
      facet_wrap(~as.factor(map)) +
      geom_line(data = marey, aes(x = .data$physicalPosition/10^6, y = .data$geneticPositioncM), colour = "black") +
      geom_ribbon(data = marey, aes(x = .data$physicalPosition/10^6, y = .data$geneticPositioncM, ymin = .data$lowerGeneticPositioncM, ymax = .data$upperGeneticPositioncM),
                  fill = "darkorange", colour = "darkorange", alpha = 0.3) +
      labs(x = "Genomic position (Mb)", y = "Genetic distance (cM)") +
      theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(color="black", face="bold.italic",hjust = 0.5),
            plot.subtitle = element_text(color="black",hjust = 0.5),
            axis.title.x = element_text(color="black"),
            axis.title.y = element_text(color="black"),
            axis.text=element_text(colour="black"),
            legend.position = "right")
  } else {
    p = ggplot2::ggplot(data = df, aes(x = .data$phys/10^6, y = .data$gen, colour = .data$vld)) +
      geom_point() +
      scale_colour_manual(values=c("TRUE" = "dodgerblue4", "FALSE" = "firebrick4")) +
      facet_grid(~as.factor(set)) +
      facet_wrap(~as.factor(map)) +
      labs(x = "Genomic position (Mb)", y = "Genetic distance (cM)") +
      theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(color="black", face="bold.italic",hjust = 0.5),
            plot.subtitle = element_text(color="black",hjust = 0.5),
            axis.title.x = element_text(color="black"),
            axis.title.y = element_text(color="black"),
            axis.text=element_text(colour="black"),
            legend.position = "right")
  }

  if (return.plot) {
    return(p)
  } else {
    return(df)
  }
}



#' Plot the recombination map
#'
#' @param x a Marey map object.
#' @param return.plot TRUE or FALSE, whether to return a plot or a data frame (default = TRUE returns the plot)
#'
#' @return a plot of the recombination map.
#' @export
#'
#' @import ggplot2
#'
plot_recombinationMap = function(x, return.plot = TRUE) {
  x$recMap$point = (x$recMap$start + x$recMap$end)/2

  x$recMap$recRate = x$recMap$recRate * 10^6
  x$recMap$upperRecRate = x$recMap$upperRecRate * 10^6
  x$recMap$lowerRecRate = x$recMap$lowerRecRate * 10^6

  p = ggplot2::ggplot(data = x$recMap, aes(x = .data$point/10^6, y = .data$recRate)) +
    geom_line() +
    geom_ribbon(aes(x = .data$point/10^6, ymin = .data$lowerRecRate, ymax = .data$upperRecRate), fill = "gray30", alpha = 0.2) +
    facet_grid(~as.factor(set)) +
    facet_wrap(~as.factor(map), scales = "free") +
    labs(x = "Genomic position (Mb)", y = "Recombination rate (cM/Mb)") +
    theme(axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", face="bold.italic",hjust = 0.5),
          plot.subtitle = element_text(color="black",hjust = 0.5),
          axis.title.x = element_text(color="black"),
          axis.title.y = element_text(color="black"),
          axis.text=element_text(colour="black"),
          legend.position = "right")
  if (return.plot) {
    return(p)
  } else {
    return(df)
  }
}

