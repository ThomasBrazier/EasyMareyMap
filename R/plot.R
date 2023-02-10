#' Plot the Marey maps
#'
#' @param x A Marey map object
#'
#' @return A plot of the Marey maps
#' @export
#'
#' @examples
plot.mareyMap = function(x) {
  p = ggplot(data = x$mareyMap, aes(x = phys/10^6, y = gen, col = vld)) +
    geom_point() +
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
  return(p)
}
