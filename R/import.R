##------------------------------------------------------------------------
#' import_clonedata
#'
#' Imports all clone data as list of data frames by simulation number
#'
#' @param filedir folder directory containing clonedata.txt
#'
#' @return clone_data - A list of data frames for each simulation
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' clone_data <- import_clonedata('./clonedata.txt')
#' }
import_clonedata <- function(filename = "./clonedata") {
  clone_df <- read.table(filename,
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F)

  clone_data <- split.data.frame(clone_df, f = clone_df$run)
  return(clone_data)
}

##------------------------------------------------------------------------
#' import_timedata
#'
#' Import time course data as list of data frames by simulation number
#'
#' @param filedir folder directory containing clonedata.txt
#'
#'
#' @return time_data - A list of data frames for each simulation
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' time_data <- import_timedata('./timedata.txt')
#' }
import_timedata <- function(filename = "./timedata") {
  time_df <- read.table(filename,
                        sep = "\t",
                        header = T,
                        stringsAsFactors = F)

  time_data <- split.data.frame(time_df, f = time_df$run)
  return(time_data)
}


##------------------------------------------------------------------------
#' Import summary of simulations and convert to a list
#'
#' @param filedir folder directory containing clonedata.txt
#'
#'
#' @return sim_data - A list of variables for all simulations
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' time_data <- import_simdata('./sim_stats.txt')
#' }
import_simdata <- function(filename = "./sim_stats.txt") {
  sim_df <- read.csv(filename,
                         header = F,
                         stringsAsFactors = F,
                         row.names = 1,
                         col.names = c("name", "value"))
  sim_data <- setNames(split(sim_df$value, seq(nrow(sim_df))), rownames(sim_df))
  return(sim_data)
}

##------------------------------------------------------------------------
#' Import single cell samples
#'
#' @param filedir folder directory containing clonedata.txt
#'
#'
#' @return sim_data - A list of variables for all simulations
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' sample_data <- import_sampledata('./sampledata.txt')
#' }
import_sampledata <- function(filename = "./sampledata.txt") {
  sample_df <- read.delim(filename,
                     header = T, sep = "\t",
                     stringsAsFactors = F)
  sample_df <- split.data.frame(sample_df, f = sample_df$run)
  return(sample_df)
}

##------------------------------------------------------------------------
#' import_siapop
#'
#' Import all data from a simulation
#'
#' The structure of the list contains the data in the first element and
#' simulation data in the second. The data element is a list of each run with
#' clone data, time data, and samples as elements in each.
#'
#' Either \code{filedir} is required and will import the files from the folder
#' or individual filenames need to be supplied. The clone data file is required
#' but if other files are not supplied they are given as NULL values.
#'
#' @param filedir folder directory containing siapop output
#' @param clonedata_file file of clone data file (required if \code{filedir} not supplied.
#' @param timedata_file file of time course data file
#' @param sampledata_file file of sample data file
#' @param simdata_file file of simulation information
#'
#'
#'
#' @return siapop_data Complete output of siapop in list format.
#' @export
#' @examples
#' \dontrun{
#' siapopConstant(seed = 17, max_pop = 1000, mutation_prob = 0.05,
#'                observation_frequency = 1, detection_threshold = 0.005,
#'                num_samples = 1, sample_size = 100)
#' siapop_data <- import_siapop('./')
#' }
import_siapop <- function(filedir, clonedata_file = NULL, timedata_file = NULL,
                          sampledata_file = NULL, simdata_file = NULL) {

  if(missing(filedir)){
    if(!is.null(clonedata_file) && file.exists(clonedata_file)){
      clone_data <- import_clonedata(clonedata_file)
    }
    else{
      stop("Need clone data file. File not supplied or does not exist")
    }

    if(!is.null(timedata_file) && file.exists(timedata_file)){
      time_data <- import_timedata(timedata_file)
      time_data <- time_data[names(time_data) %in% names(clone_data)]
    }
    else{
      time_data <- vector("list", length(clone_data))
    }

    if(!is.null(sampledata_file) && file.exists(sampledata_file)){
      sample_data <- import_sampledata(sampledata_file)
    }
    else{
      sample_data <- vector("list", length(clone_data))
    }

    if(!is.null(simdata_file) && file.exists(simdata_file)){
      sim_data <-import_simdata(simdata_file)
    }
    else{
      sim_data <-vector("list", length(clone_data))
    }
  }
  else{

    if(file.exists(paste(filedir, "/clonedata.txt", sep = ""))){
      clone_data <- import_clonedata(paste(filedir, "/clonedata.txt", sep = ""))
    }
    else{
      stop("Need clone data file. File 'clonedata.txt' does not exist at directory location")
    }
    if(file.exists(paste(filedir, "/timedata.txt", sep = ""))){
      time_data <- import_timedata(paste(filedir, "/timedata.txt", sep = ""))
      time_data <- time_data[names(time_data) %in% names(clone_data)]
    }
    else{
      time_data <- vector("list", length(clone_data))
    }
    if(file.exists(paste(filedir, "/sampledata.txt", sep = ""))){
      sample_data <- import_sampledata(paste(filedir, "/sampledata.txt", sep = ""))
    }
    else{
      sample_data <- vector("list", length(clone_data))
    }

    if(file.exists(paste(filedir, "/sim_stats.txt", sep = ""))){
      sim_data <- import_simdata(paste(filedir, "/sim_stats.txt", sep = ""))
    }
    else{
      sim_data <- vector("list", length(clone_data))
    }
  }

  siapop_data <- vector("list", length(clone_data))
  names(siapop_data) <- names(clone_data)
  for(i in 1:length(clone_data)){
    siapop_data[[i]] <- list(clone_data = clone_data[[i]],
                             time_data = time_data[[i]],
                             sample_data = sample_data[[i]])
  }
  siapop_data <- list(data = siapop_data, simulation_info = sim_data)
  return(siapop_data)
}
