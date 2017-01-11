##------------------------------------------------------------------------
#' Import clone data as list of data frames by simulation number
#'
#' @param filedir folder directory containing clonedata.txt
#'
#'
#' @return clone_data - A list of data frames for each simulation
#' @export
import_clonedata <- function(filedir) {
  clone_df <- read.table(paste(filedir, "/clonedata.txt", sep = ""),
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F)

  clone_data <- split.data.frame(clone_df, f = clone_df$run)
  return(clone_data)
}

##------------------------------------------------------------------------
#' Import time course data as list of data frames by simulation number
#'
#' @param filedir folder directory containing clonedata.txt
#'
#'
#' @return time_data - A list of data frames for each simulation
#'
import_timedata <- function(filedir) {
  time_df <- read.table(paste(filedir, "/timedata.txt", sep = ""),
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
import_simdata <- function(filedir) {
  sim_df <- read.csv(paste(filedir, "/sim_stats.txt", sep = ""),
                         header = F,
                         stringsAsFactors = F,
                         row.names = 1,
                         col.names = c("name", "value"))
  sim_data <- setNames(split(sim_df$value, seq(nrow(sim_df))), rownames(sim_df))
  return(sim_data)
}

##------------------------------------------------------------------------
#' Import all data from a simulation and outputs a list
#'
#' @param filedir folder directory containing siapop output
#'
#'
#' @return siapop_data - Complete output of siapop in list format
#' @export
import_siapop <- function(filedir) {
  clone_data <- import_clonedata(filedir)
  time_data <- import_timedata(filedir)
  sim_data <- import_simdata(filedir)

  siapop_data <- list(clone_data = clone_data, time_data = time_data,
                      simulation_stats = sim_data)
}
