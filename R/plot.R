#' Plot the Marey maps
#'
#' @param x a Marey map object.
#'
#' @return a plot of the Marey maps.
#' @export
#'
#' @import ggplot2
#'
plot_marey = function(x) {
  df = x$mareyMap
  # Add the Marey interpolated function
  marey = x$mareyCI
  if (nrow(marey) > 0) {
    marey$vld = TRUE

    p = ggplot2::ggplot(data = df, aes(x = phys/10^6, y = gen, colour = vld)) +
      geom_point(alpha = 0.4) +
      scale_colour_manual(values=c("TRUE" = "dodgerblue4", "FALSE" = "firebrick4")) +
      facet_grid(~as.factor(set)) +
      facet_wrap(~as.factor(map)) +
      geom_line(data = marey, aes(x = physicalPosition/10^6, y = geneticPositioncM), colour = "black") +
      geom_ribbon(data = marey, aes(x = physicalPosition/10^6, y = geneticPositioncM, ymin = lowerGeneticPositioncM, ymax = upperGeneticPositioncM),
                  fill = "gray30", colour = "gray60", alpha = 0.2) +
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
    p = ggplot2::ggplot(data = df, aes(x = phys/10^6, y = gen, colour = vld)) +
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

  return(p)
}
