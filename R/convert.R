##------------------------------------------------------------------------
#' create_sample_adj_matrix
#'
#' Creates an adjacency matrix from a sample of SIApop. The adjacency matrix
#' assumes rows are individuals and columns are alleles and returns a TRUE
#' in the i,j element of the matrix if individual sample i contains allele j.
#'
#' The adjacency matrix create allows conversion a distance matrix for the
#' creation of phylogenetic trees with package \code{ape} or converted to a
#' \code{phyDat} object in the \code{phanghorn} package.
#'
#' @param samp SIApop single sample data frame containing a column for the
#'  unique id (\code{unique_id}) and number of samples for each id
#'  (\code{number_obs}).
#'
#' @return create_sample_adj_matrix - returns an adjacency matrix for a sample
#' @export
#' @examples
#' \dontrun{
#' # simulate process and import
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' samp_data <- import_sampledata('./sampledata.txt')$`1`
#' create_sample_adj_matrix(samp_data)
#' }
create_sample_adj_matrix <- function(samp){
  if(length(unique(samp$sample_number)) > 1){
    warning("More than 1 sample present. Consider subsetting")
  }

  names <- rep(samp$unique_id, samp$number_obs)
  alleles <- unique(unlist(strsplit(names, ">")))

  adj_matrix <- sapply(alleles, function(allele) sapply(strsplit(names,">"),
    function(idlist)  allele %in% idlist))
  row.names(adj_matrix) <- make.names(names, unique = TRUE)
  return(adj_matrix)
}

##------------------------------------------------------------------------
#' create_adj_matrix
#'
#' Creates an adjacency matrix from a population of SIApop. The adjacency matrix
#' assumes rows are individuals and columns are alleles and returns a TRUE
#' in the i,j element of the matrix if individual sample i contains allele j.
#'
#' The adjacency matrix create allows conversion a distance matrix for the
#' creation of phylogenetic trees with package \code{ape} or converted to a
#' \code{phyDat} object in the \code{phanghorn} package.
#'
#' @param clone_data SIApop single simulation data frame containing a column
#'  for the unique id (\code{unique_id}) and number of samples for each id
#'  (\code{number_obs}).
#'
#' @return create_adj_matrix - returns an adjacency matrix for a sample
#' @export
#' @examples
#' \dontrun{
#' # simulate process and import
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' clone_data <- import_clonedata('./clonedata.txt')$`1`
#' create_adj_matrix(clone_data)
#' }
create_adj_matrix <- function(clone_data){
  names <- clone_data$unique_id
  alleles <- unique(unlist(strsplit(names, ">")))
  adj_matrix <- sapply(alleles, function(allele) sapply(strsplit(names,">"),
    function(idlist)  allele %in% idlist))
  row.names(adj_matrix) <- names
  return(adj_matrix)
}


##------------------------------------------------------------------------
#' convert_fishplot
#'
#' Converts time data from a data frame to a fish object in package 'fishplot'
#' for plotting.
#'
#' @param time_data dataframe of time course data
#' @param timepoints timepoints to use
#' @param threshold minimum allele frequency to count a clone at
#'
#' @return convert_fishplot - returns a fish object for use with fishPlot()
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' time_data <- import_timedata('./timedata.txt')$`1`
#' fish <- convert_fishplot(time_data, timepoints = 1:10)
#' fishPlot(fish)
#' }
convert_fishplot <- function(time_data, threshold = 0.0001, timepoints = NULL){

  if (requireNamespace("fishplot", quietly = FALSE)) {

    if(is.null(timepoints)) timepoints <- unique(time_data$time)

    max_cell_count <- max((time_data %>% group_by(time) %>%
                           summarize(numcells = sum(numcells)))$numcells)

    time_data <- time_data %>% filter(time %in% timepoints) %>%
      select(time, unique_id, allelefreq)


    to_keep <- (time_data %>% mutate(allelefreq =
                                              allelefreq / max_cell_count) %>%
                  group_by(unique_id) %>%
                  summarize(maxfreq = max(allelefreq)) %>%
                  filter(maxfreq > threshold))$unique_id

    parents <- .match_parents(to_keep)

    time_data <- time_data %>% filter(unique_id %in% to_keep) %>%
      mutate(allelefreq = allelefreq / max_cell_count * 100)

    frac.table <- as.matrix((time_data %>%
                               tidyr::spread(time, allelefreq, fill = 0))[,-1])

    fish = fishplot::createFishObject(frac.table, parents,
                                      timepoints = timepoints,
                                      col = rainbow(nrow(frac.table)),
                                      clone.labels = to_keep)
    fish = fishplot::layoutClones(fish)
    return(fish)
  }
  else{
    stop("package 'fishplot' can be installed with
          devtools::install_github(\"chrisamiller/fishplot\")")
  }
}

##------------------------------------------------------------------------
#' convert_ggmuller
#'
#' Converts time data to a Muller_df object in package 'ggmuller' for plotting.
#'
#' @param time_data dataframe of time course data
#' @param timepoints timepoints to use
#' @param threshold minimum allele frequency to count a clone at
#' @param freqplot if TRUE then plot as a frequency relative to each time point
#' @param reduce if TRUE then reduce clone names to their final ID rather than
#'  their ancestry
#'
#'
#' @return convert_ggmuller - returns a Muller_df tibble for Muller_plot
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#'  time_data <- import_timedata('./timedata.txt')$`1`
#'  mullerdf <- convert_ggmuller(time_data)
#'  ggmuller::Muller_plot(mullerdf)
#' }
convert_ggmuller <- function(time_data, threshold = 0.001, timepoints = NULL,
                             freqplot = FALSE, reduce = TRUE, ...){

  if(is.null(timepoints)) timepoints <- unique(time_data$time)
  time_data <- time_data %>% filter(time %in% timepoints)

  ##################################################################
  ### if want to call allelefreq as freq w.r.t total number of cells
  ### throughout time instead of at a single time point
  # max_cell_count <- max((time_data %>% group_by(time) %>%
  #                          summarize(numcells = sum(numcells)))$numcells)
  # to_keep <- (time_data %>%
  #               mutate(allelefreq = allelefreq / max_cell_count) %>%
  #               group_by(unique_id) %>%
  #               summarize(maxfreq = max(allelefreq))%>%
  #               filter(maxfreq >= threshold))$unique_id
  ##################################################################

  to_keep <- (time_data %>% group_by(time) %>%
                mutate(allelefreq = allelefreq / sum(numcells)) %>%
                ungroup %>% group_by(unique_id) %>%
                summarize(maxfreq = max(allelefreq)) %>%
                filter(maxfreq >= threshold))$unique_id


  edgelist <- .create_edge.list(to_keep, reduce = reduce)

  pop_df <- time_data %>% filter(unique_id %in% to_keep) %>%
    rename(Generation = time, Identity = unique_id,
                  Population = numcells)

  if(reduce) pop_df$Identity <- sapply(pop_df$Identity, .pop, ">")

  max_cell_count <- max((pop_df %>% group_by(Generation) %>%
                           summarize(numcells =
                                              sum(Population)))$numcells)


  dummy_ancestor <- data.frame(Generation = -1, Identity = "0",
                               Population = 0, stringsAsFactors = F)
  pop_df <- bind_rows(dummy_ancestor, pop_df)
  pop2 <- pop_df
  pop_df <- pop_df %>% select(Generation, Identity, Population)

  pop_df <- pop_df %>% tidyr::spread(Identity, Population, fill = 0)

  if(freqplot) {
    pop_df$'0' <- 0
  } else {
    pop_df$'0' <- max_cell_count - rowSums(pop_df[,-1]) + 1e-8
  }

  pop_df <- pop_df %>% tidyr::gather(Identity, "Population", -1) %>%
    arrange(Generation)
  pop_df <- pop_df %>% left_join(pop2)

  if (requireNamespace("ggmuller", quietly = FALSE)) {
    return(ggmuller::get_Muller_df(edgelist, pop_df))
  }
  else{
    stop("package 'ggmuller' can be installed with
         devtools::install_github(\"robjohnnoble/ggmuller\")")
  }
}


##------------------------------------------------------------------------
#' convert_igraph
#'
#' Converts data to igraph in tree layout to view the gene tree at the final
#' time or a specific time (if importing the time data)
#'
#' @param clonedata data frame of a single timepoint
#' @param threshold minimum allele frequency to count a clone at
#' @param size "count" makes the size of nodes proportional to clone sizes.
#' @param color = c("fitness", "age", "count") the color of nodes
#'
#' @return convert_igraph - returns an igraph object
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' clone_data <- import_clonedata('./clonedata.txt')$`1`
#' clonegraph <- convert_igraph(clone_data)
#' plot(clonegraph)
#' }
convert_igraph <- function(clonedata, threshold = 0.01, size = NULL, color = NULL){
  if (requireNamespace("igraph", quietly = FALSE)) {

    to_keep <- (clonedata %>%
                  mutate(allelefreq = allelefreq / sum(allelefreq)) %>%
                  filter(allelefreq >= threshold))$unique_id

    clonedata <- clonedata %>% filter(unique_id %in% to_keep)

    edgelist <- .create_edge.list(to_keep)
    edgelist <- edgelist[edgelist$Parent != 0,]

    birth_times <- sapply(edgelist$Identity, .parent_age_at_split, edgelist,
                          clonedata$initialtime)
    nodes <- sapply(clonedata$unique_id, .pop, ">")
    birth_times <- max(clonedata$initialtime) - clonedata$initialtime

    genetree_graph <- igraph::graph.data.frame(edgelist)
    genetree_graph <- igraph::add_layout_(genetree_graph, igraph::as_tree())
    igraph::graph_attr(genetree_graph, "layout")[,2] <-
      birth_times[match(igraph::as_ids(igraph::V(genetree_graph)), nodes)]

    if(!is.null(color)){
      if(color == "fitness"){
        value <- clonedata$birthrate - clonedata$deathrate
        value <- value[match(igraph::as_ids(igraph::V(genetree_graph)), nodes)]
      } else if(color == "age"){
        value <- clonedata$initialtime
        value <- value[match(igraph::as_ids(igraph::V(genetree_graph)), nodes)]

      } else if(color == "count"){
        value <- clonedata$numcells
        value <- value[match(igraph::as_ids(igraph::V(genetree_graph)), nodes)]
      } else{
        stop("no applicable variable for color")
      }

      num_sequence <- seq(min(value), max(value), length.out = 100)
      col_sequence <- colorRampPalette(c("blue", "red"))(100)
      cols <- sapply(value, function(x) which.min(abs(num_sequence - x)))

      igraph::vertex_attr(genetree_graph, "color") <- col_sequence[cols]
    }

    if(!is.null(size)){
      if(size == "count"){
        value <- log10(clonedata$numcells)
        value <- value[match(igraph::as_ids(igraph::V(genetree_graph)), nodes)]
      } else {
        stop("no applicable variable for size")
      }
      num_sequence <- seq(min(value), max(value), length.out = 5)
      sizes <- sapply(value, function(x) which.min(abs(num_sequence - x)))
      igraph::vertex_attr(genetree_graph, "size") <- num_sequence[sizes]
    }
    return(genetree_graph)
  }
}
