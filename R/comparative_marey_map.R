#' A meta Marey map object for comparative analyses
#' 
#' This function instantiate a meta marey map object storing mutliple `marey_map` objects
#' A data frame must be provided with one or many datasets and one or many chromosomes
#' Marey map functions (e.g. metrics, recombination map estimaion) can be applied in batch
#' on this object to perform comparative analyses
#'
#' @param x a data frame of a Marey map with mandatory columns "set",	"map",	"mkr",	"phys",	"gen".	"vld" is optional and added if not provided.
#' @param chromosome_length (optional) a vector of chromosome length. If not provided, max genomic position is assumed to be chromosome end
#' 
#' 
#' @slot data a list of `marey_map` objects
#' @slot set a vector of dataset names
#' @slot map a vector of chromosome names
#' 
#' @return a 'comparative_marey_map' object.
#' @export
#'
comparative_marey_map = function(x = data.frame(),
                                chromosome_length = numeric()) {
  x$set = as.factor(x$set)
  x$map = as.factor(x$map)
  x$mkr = as.character(x$mkr)
  x$phys = as.numeric(x$phys)
  x$gen = as.numeric(x$gen)
  if (!('vld' %in% colnames(x))) {
    x$vld = logical()
  }
  
  l = split(x, x[,c('set', 'map')])
  set = unlist(lapply(l, function(x) as.character(unique(x$set))))
  map = unlist(lapply(l, function(x) as.character(unique(x$map))))
  
  l = lapply(l, function(x) marey_map(x, chromosome_length = chromosome_length))
  
  df = list(data = l,
            set = set,
            map = map)
  
  structure(df,
            class = "comparative_marey_map"
      )
}




#' Estimate recombination maps in meta Marey map object
#'
#' @param x a `comparative_marey_map` object
#' @param method the interpolation method to apply
#' @param verbose whether to print messages and progress bar
#' @param ... additional arguments
#' 
#' @return a `comparative_marey_map` object with recombination maps updated
#' @export
#'
comparative_recombination_maps = function(x,
                                          method = 'loess',
                                          verbose = TRUE, ...) {
  l = x$data
  
  x$data = lapply(l, function(x) recombination_map(x, method = method, verbose = verbose))
  
  return(x)
}



#' Compute summary statistics on a meta `comparative_marey_map` object
#'
#' @param x a `comparative_marey_map` object
#' @param statistics a vector of statistics to compute
#' @param ... additional arguments
#' 
#' @slot set a vector of dataset names
#' @slot map a vector of chromosome names
#' @slot statistics one slot per statistic computed, see details for the complete list of statistics
#' 
#' @details
#' The complete list of statistics that can be computed is
#' `mean`, `median`, `weighted.mean`, `variance`, `gini`, `peripherybias` and `coefficient_variation`
#' 
#' @importFrom stats aggregate
#' @importFrom stats var
#' 
#' @return a list of summary statistics
#' @export
#'
compute_stats_marey = function(x,
                               statistics = c('mean', 'median'), ...) {
  list_stats = list()
  list_stats$set = as.character(x$set)
  list_stats$map = as.character(x$map)
  
  # marey = comparative_marey_to_dataframe(x)
  recmap = comparative_recmap_to_dataframe(x)
  
  if ('mean' %in% statistics) {
    m = as.numeric(aggregate(recRate ~ set + map, recmap, mean)$recRate)
    list_stats$mean = m
  }
  if ('median' %in% statistics) {
    m = as.numeric(aggregate(recRate ~ set + map, recmap, median)$recRate)
    list_stats$median = m
  }
  if ('weighted.mean' %in% statistics) {
    wm = c()
    for (i in 1:length(x$data)) {
      df = x$data[i]
      res = weighted.mean(df[[1]])
      wm = c(wm, res)
    }
    list_stats$weighted.mean = wm
  }
  if ('variance' %in% statistics) {
    v = as.numeric(aggregate(recRate ~ set + map, recmap, var)$recRate)
    list_stats$variance = v
  }
  if ('gini' %in% statistics) {
    g = c()
    for (i in 1:length(x$data)) {
      df = x$data[i]
      res =gini(df[[1]])
      g = c(g, res)
    }
    list_stats$gini = g
  }
  if ('peripherybias' %in% statistics) {
    pbr = c()
    for (i in 1:length(x$data)) {
      df = x$data[i]
      res = peripherybias(df[[1]])
      pbr = c(pbr, res)
    }
    list_stats$peripherybias = pbr
  }
  if ('coefficient_variation' %in% statistics) {
    cv = c()
    for (i in 1:length(x$data)) {
      df = x$data[i]
      res = coefficient_variation(df[[1]])
      cv = c(cv, res)
    }
    list_stats$coefficient_variation = cv
  }

  return(list_stats)
}




#' Plot comparative Marey maps
#'
#' @param x a `comparative_marey_map` object
#' @param group the grouping factor in `ggplot`, either `set`, `map` or `set + map`
#' @param relative_distance_x logical, whether to plot the X axis as relative distances (0,1) to compare different genomes
#' @param relative_distance_y logical, whether to plot the Y axis as relative distances (0,1) to compare different genetic maps
#' @param ... arguments passed to the generic summary function.
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @importFrom stats as.formula
#' 
#' @return a `ggplot2` object of Marey maps
#' @export
#'
plot_comparative_marey = function(x,
                                  group = 'set + map',
                                  relative_distance_x = FALSE,
                                  relative_distance_y = FALSE, ...) {
  
  df = comparative_marey_to_dataframe(x)
  
  scaling = 10^6 # Convert bp to Mb
  
  # Relative distances
  if (relative_distance_x) {
    df$max_phys = unlist(lapply(1:nrow(df), function(x) {max(df$phys[which(df$set == df$set[x] & df$map == df$map[x])], na.rm = T)}))
    df$phys = df$phys / df$max_phys
    scaling = 1
    lab = labs(x = "Relative genomic position", y = "Genetic distance (cM)")
  }
  if (relative_distance_y) {
    df$max_gen = unlist(lapply(1:nrow(df), function(x) {max(df$gen[which(df$set == df$set[x] & df$map == df$map[x])], na.rm = T)}))
    df$gen = df$gen / df$max_gen
    lab = labs(x = "Genomic position (Mb)", y = "Relative genetic distance")
  }
  if (relative_distance_x & relative_distance_y){
    lab = labs(x = "Relative genomic position", y = "Relative genetic distance")
  }
  if (!relative_distance_x & !relative_distance_y){
    lab = labs(x = "Genomic position (Mb)", y = "Genetic distance (cM)")
  }

  # Add the Marey interpolated function
  marey = comparative_interpolation_to_dataframe(x)
  if (relative_distance_x) {
    marey$max_phys = unlist(lapply(1:nrow(marey), function(x) {max(marey$physicalPosition[which(marey$set == marey$set[x] & marey$map == marey$map[x])], na.rm = T)}))
    marey$physicalPosition = marey$physicalPosition / marey$max_phys
  }
  if (relative_distance_y) {
    marey$max_gen = unlist(lapply(1:nrow(marey), function(x) {max(marey$upperGeneticPositioncM[which(marey$set == marey$set[x] & marey$map == marey$map[x])], na.rm = T)}))
    marey$geneticPositioncM = marey$geneticPositioncM / marey$max_gen
    marey$upperGeneticPositioncM = marey$upperGeneticPositioncM / marey$max_gen
    marey$lowerGeneticPositioncM = marey$lowerGeneticPositioncM / marey$max_gen 
  }
  
  
  if (group == 'set') {
    grouping = as.formula(~as.factor(.data$set))
    facet = facet_wrap(grouping, scales = "free")
    point_rec = geom_point(aes(colour = as.factor(.data$map), fill = as.factor(.data$map)), alpha = 0.2)
    line_rec = geom_line(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM, group = as.factor(.data$map)), fill = "black")
    ribbon_rec = geom_ribbon(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM, ymin = .data$lowerGeneticPositioncM, ymax = .data$upperGeneticPositioncM, group = as.factor(.data$map)),
                             alpha = 0.4)

  }
  if (group == 'map') {
    grouping = as.formula(~ as.factor(.data$map))
    facet = facet_wrap(grouping, scales = "free")
    point_rec = geom_point(aes(colour = as.factor(.data$set), fill = as.factor(.data$set)), alpha = 0.2)
    line_rec = geom_line(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM, group = as.factor(.data$set)), fill = "black")
    ribbon_rec = geom_ribbon(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM, ymin = .data$lowerGeneticPositioncM, ymax = .data$upperGeneticPositioncM, group = as.factor(.data$set)),
                             alpha = 0.4)
  }
  if (group == 'set + map') {
    grouping = as.formula(~as.factor(.data$map) + as.factor(.data$set))
    facet = facet_grid(grouping, scales = "free")
    point_rec = geom_point(alpha = 0.2)
    line_rec = geom_line(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM), colour = "black")
    ribbon_rec = geom_ribbon(data = marey, aes(x = .data$physicalPosition/scaling, y = .data$geneticPositioncM, ymin = .data$lowerGeneticPositioncM, ymax = .data$upperGeneticPositioncM),
                               fill = "darkorange", colour = "darkorange", alpha = 0.3)
  }
  
  
  if (nrow(marey) > 0) {
    marey$vld = TRUE
    
    p = ggplot2::ggplot(data = df, aes(x = .data$phys/scaling, y = .data$gen)) +
      point_rec +
      line_rec +
      ribbon_rec +
      facet +
      lab +
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
    p = ggplot2::ggplot(data = df, aes(x = .data$phys/scaling, y = .data$gen)) +
      point_rec +
      facet +
      labs(x = "Genomic position (Mb)", y = "Genetic distance (cM)", colour = "dataset") +
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
  
  p
}



#' Plot comparative recombination maps
#'
#' @param x a `comparative_marey_map` object
#' @param group the grouping factor in `ggplot`, either `set`, `map` or `set + map`
#' @param relative_distance_x logical, whether to plot the X axis as relative distances (0,1) to compare different genomes
#' @param relative_distance_y logical, whether to plot the Y axis as relative distances (0,1) to compare maps with different average recombination rates
#' 
#' @return a `ggplot2` object of Marey maps
#' @export
#' 
#' @import ggplot2
#'
plot_comparative_recmap = function(x,
                                   group = 'set + map',
                                   relative_distance_x = FALSE,
                                   relative_distance_y = FALSE) {
  
  x = comparative_recmap_to_dataframe(x)

  x$point = (x$start + x$end)/2
  
  x$recRate = x$recRate * 10^6
  x$upperRecRate = x$upperRecRate * 10^6
  x$lowerRecRate = x$lowerRecRate * 10^6
  
  # Relative distances
  if (relative_distance_x) {
    x$max_phys = unlist(lapply(1:nrow(x), function(y) {max(x$end[which(x$set == x$set[y] & x$map == x$map[y])], na.rm = T)}))
    x$point = x$point / x$max_phys
    lab = labs(x = "Relative genomic position", y = "Recombination rate (cM/Mb)")
  }
  if (relative_distance_y) {
    x$max_rec = unlist(lapply(1:nrow(x), function(y) {max(x$upperRecRate[which(x$set == x$set[y] & x$map == x$map[y])], na.rm = T)}))
    x$recRate = x$recRate / x$max_rec
    x$upperRecRate = x$upperRecRate / x$max_rec
    x$lowerRecRate = x$lowerRecRate / x$max_rec
    lab = labs(x = "Genomic position (Mb)", y = "Relative recombination rate")
  }
  if (relative_distance_x & relative_distance_y){
    lab = labs(x = "Relative genomic position", y = "Recombination rate")
  }
  if (!relative_distance_x & !relative_distance_y){
    lab = labs(x = "Genomic position (Mb)", y = "Recombination rate (cM/Mb)")
  }
  
  if (group == 'set') {
    grouping = as.formula(~as.factor(.data$set))
    facet = facet_wrap(grouping, scales = "free")
    line_rec = geom_line(aes(group = .data$map, color = .data$map))
    ribbon_rec = geom_ribbon(aes(x = .data$point/10^6, ymin = .data$lowerRecRate, ymax = .data$upperRecRate, fill = .data$map), alpha = 0.2)
  }
  if (group == 'map') {
    grouping = as.formula(~as.factor(.data$map))
    facet = facet_wrap(grouping, scales = "free")
    line_rec = geom_line(aes(group = .data$set, color = .data$set))
    ribbon_rec = geom_ribbon(aes(x = .data$point/10^6, ymin = .data$lowerRecRate, ymax = .data$upperRecRate, fill = .data$set), alpha = 0.2)
  }
  if (group == 'set + map') {
    grouping = as.formula(~as.factor(.data$map) + as.factor(.data$set))
    facet = facet_grid(grouping, scales = "free")
    line_rec = geom_line()
    ribbon_rec = geom_ribbon(aes(x = .data$point/10^6, ymin = .data$lowerRecRate, ymax = .data$upperRecRate), fill = 'darkgray', alpha = 0.2)
  }
  
  p = ggplot2::ggplot(data = x, aes(x = .data$point/10^6, y = .data$recRate)) +
    line_rec +
    ribbon_rec +
    # facet_grid(~as.factor(set)) +
    # facet_wrap(grouping, scales = "free") +
    facet +
    lab +
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
  
  p
}



#' Merge `comparative_marey_map` objects
#'
#' @param x a list of `comparative_marey_map` objects to be merged
#' 
#' @return a `comparative_marey_map` object of Marey maps
#' @export
#'
merge_comparative_marey = function(x = list()) {
  
  new_object = comparative_marey_map()
  
  n_list = length(x)
  
  for (i in 1:n_list) {
    d = x[[i]]
    new_object$data = c(new_object$data, d$data)
    new_object$set = c(new_object$set, d$set)
    new_object$map = c(new_object$map, d$map)
  }
  
  # check consistency
  stopifnot(length(new_object$set) == length(new_object$data))
  stopifnot(length(new_object$set) == length(new_object$map))
  
  return(new_object)
}


#' Subset a `comparative_marey_map` object
#'
#' @param x a `comparative_marey_map` object
#' @param subset logical expression indicating `set` and `map` to keep
#' 
#' @return a `comparative_marey_map` object of Marey maps
#' @export
#'
subset_comparative_marey = function(x,
                                    subset = c()) {
  
  x2 = comparative_marey_map()
  x2$data = subset(x$data, subset)
  x2$set = subset(x$set, subset)
  x2$map = subset(x$map, subset)
  
  return(x2)
}




#' Summary of a `comparative_marey_map` object
#'
#' @param object a `comparative_marey_map` object to summarize
#' @param ... additional arguments
#' 
#' @return a summary
#'
#' @method summary comparative_marey_map
#' @export
#' 
summary.comparative_marey_map = function(object, ...) {
  dataset = paste0(c(as.character(head(object$set, 3)), "..."))
  n_maps = length(object$set)
  length_linkage_map = unlist(lapply(object$data, function(x) max(x[[1]]$gen)))
  length_genome_Mb = unlist(lapply(object$data, function(x) max(x[[1]]$phys)))
  length_genome_Mb = length_genome_Mb / 1000000
  n_markers = unlist(lapply(object$data, function(x) length(x[[1]]$gen)))
  density_markers = n_markers / length_genome_Mb
  
  cat("============== Summary of the comparative marey map ==============\n",
      "Datasets: ", dataset, "\n",
      "Number of maps: ", n_maps, "\n",
      "Mean linkage map length (cM): ", mean(length_linkage_map, na.rm = TRUE), "\n",
      "Mean chromosome size (Mb): ", mean(length_genome_Mb, na.rm = TRUE), "\n",
      "Mean number of markers: ", mean(n_markers, na.rm = TRUE), "\n",
      "Mean marker density (marker/Mb): ", mean(density_markers, na.rm = TRUE), "\n",
      "==================================================================\n")
}






#' Convert Marey maps in a `comparative_marey_map` objects to `data.frame`
#'
#' @param x a `comparative_marey_map` object
#' 
#' @return a `data.frame` with `marey_map` Marey maps concatenated
#' @export
#'
comparative_marey_to_dataframe = function(x) {
  list_df = list()
  
  for (i in 1:length(x$data)) {
    df = x$data[i]
    
    out_df = df[[1]]$mareyMap

    out_df = list(out_df)
    list_df = append(list_df, out_df)
  }
  df = do.call("rbind", list_df)

  return(df)
}


#' Convert Marey interpolated functions in a `comparative_marey_map` objects to `data.frame`
#'
#' @param x a `comparative_marey_map` object
#' 
#' @return a `data.frame` with `marey_map` interpolation concatenated
#' @export
#'
comparative_interpolation_to_dataframe = function(x) {
  list_df = list()
  
  for (i in 1:length(x$data)) {
    df = x$data[i]
    
    out_df = df[[1]]$mareyCI

    out_df = list(out_df)
    list_df = append(list_df, out_df)
  }
  df = do.call("rbind", list_df)
  
  return(df)
}


#' Convert recombination maps in a `comparative_marey_map` objects to `data.frame`
#'
#' @param x a `comparative_marey_map` object
#' 
#' @return a `data.frame` with `marey_map` recombination maps concatenated
#' @export
#'
comparative_recmap_to_dataframe = function(x) {
  list_df = list()
  
  for (i in 1:length(x$data)) {
    df = x$data[i]
    df = df[[1]]$recMap
    df = list(df)
    list_df = append(list_df, df)
  }
  df = do.call("rbind", list_df)
  
  return(df)
}
