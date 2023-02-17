#' Plot the Marey maps
#'
#' @param x A Marey map object
#'
#' @return A plot of the Marey maps
#' @export
#'
#' @import ggplot2
#'
#' @examples recombinationPlot(x)
recombinationPlot = function(x) {
  x$recMap$point = (x$recMap$start + x$recMap$end)/2

  x$recMap$recRate = x$recMap$recRate * 10^6
  x$recMap$upperRecRate = x$recMap$upperRecRate * 10^6
  x$recMap$lowerRecRate = x$recMap$lowerRecRate * 10^6

  p = ggplot2::ggplot(data = x$recMap, aes(x = point/10^6, y = recRate)) +
    geom_line() +
    geom_ribbon(aes(x = point/10^6, ymin = lowerRecRate, ymax = upperRecRate), fill = "gray30", alpha = 0.2) +
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
  return(p)
}
