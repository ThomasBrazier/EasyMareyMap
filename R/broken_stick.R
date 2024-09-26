#' Plot the broken stick
#'
#' @description
#' Plot a broken stick to compare among a large set of species
#' the heterogeneity in the distribution of recombination rates along the genome
#' For a description of the method, see Brazier, T., & Glémin, S. (2022).
#' Diversity and determinants of recombination landscapes in flowering plants.
#' PLOS Genetics, 18(8), Article 8. https://doi.org/10.1371/journal.pgen.1010141
#'
#' @param x a `comparative_marey_map` object
#' @param k the number of segments (default = 10)
#' @param method the method to infer the breakpoint of segments, either `strict` to cut a marker position, or `segmented` to interpolate segments breakpoints
#' @param plot (logical) whether to plot the broken stick directly or return a data frame (default = `TRUE` will plot the figure)
#' 
#' 
#' @details
#' Broken Stick model (see White & Hill 2020 for details)
#' K segments of equal genomic size (Mb) and evaluate the genetic relative size
#' 
#' @returns a `ggplot` of the broken stick or a data frame, depending on the value of `plot`
#' 
#' @import segmented
#' @import ggplot2
#' @import reshape2
#' 
#' @export
#' 
brokenstick = function(x, k = 10, method = "strict", plot = TRUE) {
  
  # the list of set and names to process
  s = x$set
  n = x$map
  
  list_bs = list()
  
  for (i in 1:length(s)) {
    idx = (s == s[i] & n == n[i])
    # cat(s[i], n[i], "\n")
    
    subs = subset_comparative_marey(x, subset = idx)
    subs = comparative_marey_to_dataframe(subs)
    
    bs = list(brokenstick_one_map(subs, k = k, method = method))
    
    list_bs = c(list_bs, bs)
  }
  
  df = data.frame(set = s,
                  name = n)
  
  res = as.data.frame(do.call("rbind", list_bs))
  df = cbind(df, res)
  
  # Format columns in proper way
  colnames(df) = c('set', 'name', as.character(c(1:k)))
  
  # Tidy data frame
  # brokenstick = brokenstick[!(is.na(df$p1) | is.na(df$p2) | is.na(df$p3) | is.na(df$p4) | is.na(df$p5) | df(brokenstick$p6) | is.na(df$p7) | is.na(df$p8) | is.na(df$p9) | is.na(df$p10)),]
  brokenstick = melt(df)
  brokenstick$sample = paste(brokenstick$set, brokenstick$chromosome, sep = "_")
  colnames(brokenstick)=c("set","chromosome","segment","proportion.length","sample")
  
  # Set a vector of gradient color
  # Value of color is simply expected - p, the departure from the expected proportion 1/k
  brokenstick$color = (1/k) - brokenstick$proportion.length

  # Besides, estimates the ratio expected/observed (longer than expected will have lower relative recombination rate)
  brokenstick$ratio = brokenstick$proportion.length/(1/k)
  
  p = ggplot(data = brokenstick, aes(x=.data$sample, y=.data$proportion.length, fill = log10(.data$ratio)))+
    geom_bar(stat='identity', width = 1) +
    # scale_fill_manual(values = color) +
    scale_fill_viridis_c(breaks = c(-1, 0, 1), labels = c("-1", "0", "1"), direction = -1, limits = c(-1, 1),
                         values = c(0,0.7,0.8,1), option = "D") +
    # scale_fill_gradient2() +
    facet_grid(~set, scales = "free", space="free_x") +
    labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
    theme(axis.line = element_blank(),
          # axis.line.x = element_blank(), # No x axis
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", face="bold.italic", hjust = 0.5),
          plot.subtitle = element_text(color="black", hjust = 0.5),
          axis.title.x = element_text(color="black"),
          axis.title.y = element_text(color="black"),
          axis.text=element_text(colour="black"),
          axis.text.x=element_blank(), # No samples names
          axis.ticks.x=element_blank(), # No x axis
          strip.text=element_text(colour="black", angle = 90),
          legend.key = element_rect(fill = "white", linewidth = 1),
          legend.text=element_text(),
          legend.title=element_text())
  p
  
  if (plot) {
    p
  } else {
    return(brokenstick)
  }
}




#' Compute the broken stick model for a single Marey map
#' 
#' @description
#' Estimate the proportions of a broken stick model for a `comparative_marey_map` object
#' i.e. proportion of relative genetic length (cM) in k segments of equal genomic size (bp) along the chromosome
#' 
#' @param marey a single `marey_map` object
#' @param k the number of segments (default = 10)
#' @param method the method used to infer segments breakpoints
#' 
#' @importFrom stats lm
#' 
#' 
brokenstick_one_map = function(marey, k = 10, method = "strict") {
  
  marey$map = as.character(marey$map)
  marey$mkr = as.character(marey$mkr)
  marey$phys = as.numeric(as.character(marey$phys))
  
  # Create empty dataset
  stick_proportions = as.data.frame(matrix(NA, nrow = 1, ncol = k))
  colnames(stick_proportions) = paste("p", c(1:k), sep ="")
  
  if (method == "strict") {
    # Segment in k segments of equal genomic size (Mb)
    # Size of segment is total genomic size divided by k
    segment_size = round(max(marey$phys, na.rm = TRUE)/k)
    
    # A vector of segments genetic size pi (proportion of total length)
    p = numeric(k)
    for (i in 1:k) {
      # print(max(which(marey$gen < segment_size*i)))
      # p[i] = marey$phys[max(which(marey$gen < segment_size*i))]
      p[i] = max(marey$gen[which(marey$phys < segment_size*i)], na.rm = TRUE)
      # p[i] = max(marey$gen[marey$gen < segment_size*i])
      # p[i] = marey$phys[max(which(marey$gen < quantile(marey$gen, i/k)))]
      # p[i] = quantile(marey$gen, i/k)
      
      if (i > 1) {
        p[i] = p[i] - sum(p[i-1:(i-1)])
      }
    }
    # END OF METHOD STRICT
  }
  if (method == "segmented") {
    res = c()
    lin.mod <- lm(phys~gen, data = marey)
    # Estimated breakpoints in genomic distances
    res = try(segmented.mod <- segmented(lin.mod, seg.Z = ~phys, npsi=(k-1), control = seg.control(n.boot = 100, fix.npsi=TRUE)))
    # In cases of segmentation failure
    if (class(res)[1] == "try-error") {
      # Failure because not enough data
      breakpoints = c(NA, NA, NA)
    } else {
      breakpoints = c(summary(segmented.mod)$psi[,2], max(marey$phys, na.rm = TRUE))
    }
    
    for (i in 1:k) {
      p[i] = max(marey$gen[marey$phys < breakpoints[i]])
      
      if (i > 1) {
        p[i] = p[i] - sum(p[i-1:(i-1)])
      }
    }
    # END OF METHOD SEGMENTED
  }
  
  # stick_proportions is a proportion of total physical length
  stick_proportions = p/sum(p)

  return(stick_proportions)
}
