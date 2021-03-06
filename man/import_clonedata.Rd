% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{import_clonedata}
\alias{import_clonedata}
\title{import_clonedata}
\usage{
import_clonedata(filename = "./clonedata.txt")
}
\arguments{
\item{filedir}{folder directory containing clonedata.txt}
}
\value{
clone_data - A list of data frames for each simulation
}
\description{
Imports all clone data as list of data frames by simulation number
}
\examples{
\dontrun{
siapop(seed = 17, max_pop = 1000, mutation_prob = 0.05,
               observation_frequency = 1, detection_threshold = 0.005,
               num_samples = 1, sample_size = 100)
clone_data <- import_clonedata('./clonedata.txt')
}
}
