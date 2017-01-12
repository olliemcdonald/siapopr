#' collect a sample of alleles from the entire population
#'
#' @param clone_data - clone data frame
#' @param n_alleles - number of alleles to sample total
#' @param threshold - minimum frequency to allow for sampling
#'
#' @return data frame of a sample of alleles
#' @export
sample_alleles <- function(clone_data, n_alleles, threshold = 0.01){
  clone_data <- clone_data %>% select(unique_id, numcells, allelefreq)
  clone_data <- clone_data %>% mutate(alleleprop = allelefreq / sum(numcells)) %>%
    filter(alleleprop >= threshold)
  allele <- sapply(clone_data$unique_id, .pop)
  allele_sample <- rmultinom(1, n_alleles, clone_data$allelefreq)
  return(data.frame(allele = allele, reads = allele_sample))
}

#' sample cells from clone data frame (for small number of clones)
#'
#' @param clone_data - clone data frame
#' @param n_cells - number of alleles to sample total
#'
#' @return data frame containing a list of clones with respective cells sampled
#' @export
sample_cells <- function(clone_data, n_cells){
  clone_data <- clone_data %>% select(unique_id, numcells, allelefreq)
  cell_sample <- sample(rep(clone_data$unique_id, clone_data$numcells), n_cells)
  cell_sample <- data.frame(table(cell_sample), stringsAsFactors = FALSE)
  names(cell_sample) <- c("unique_id", "number_sampled")
  cell_sample$unique_id <- as.character(cell_sample$unique_id)
  clone_data <- clone_data %>% right_join(cell_sample)
  return(clone_data)
}
