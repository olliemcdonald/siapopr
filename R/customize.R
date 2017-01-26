##------------------------------------------------------------------------
#' create_fitness_template
#'
#' Creates a fitness distribution ".cpp" template file at the given location
#' which can be modified before running `compile_custom_fitness`.
#'
#' @param cppfile file name of .cpp file to create
#'
#' @export
#' @examples
#' \dontrun{
#' create_fitness_template(cppfile = "custom_dist_plugin.cpp")
#' }
create_fitness_template <- function(cppfile = "custom_dist_plugin.cpp"){
  cppfile <- unlist(strsplit(cppfile, ".", fixed = T))
  if(cppfile[length(cppfile)] != "cpp") cppfile <- c(cppfile, "cpp")
  cppfile <- paste(cppfile, collapse = ".")
  template_location <- paste(.libPaths(), "/siapopr/extras/customdist_template.cpp", sep = "")
  file.copy(template_location, cppfile)
}

##------------------------------------------------------------------------
#' create_newclone_template
#'
#' Creates a new clone plugin ".cpp" template file at the given location
#' which can be modified before running `compile_custom_newclone`.
#'
#' @param cppfile file name of .cpp file to create
#'
#' @export
#' @examples
#' \dontrun{
#' create_newclone_template(cppfile = "custom_newclone_plugin.cpp")
#' }
create_newclone_template <- function(cppfile = "custom_newclone_plugin.cpp"){
  cppfile <- unlist(strsplit(cppfile, ".", fixed = T))
  if(cppfile[length(cppfile)] != "cpp") cppfile <- c(cppfile, "cpp")
  cppfile <- paste(cppfile, collapse = ".")
  template_location <- paste(.libPaths(), "/siapopr/extras/custom_newclone_template.cpp", sep = "")
  file.copy(template_location, cppfile)
}


##------------------------------------------------------------------------
#' custom_fitness_distribution
#'
#' Allows the user to create a custom fitness distribution in C++ for use in
#' SIApop. Creates the appropriate header file and compiles the source file then
#' links into a shared object. The created shared object is used as a plugin.
#'
#' @param cppfile a character string giving the path name of a cpp file
#'
#' @export
#' @examples
#' \dontrun{
#' compile_custom_fitness(cppfile = "./plugin.cpp")
#' }
compile_custom_fitness <- function(cppfile){
  cpproot <- .pop_off(cppfile, ".", fixed = T)
  cppsuffix <- .pop(cppfile, ".", fixed = T)
  con <- file(paste(cpproot, ".h", sep = ""))

  header1 <- '#ifndef CTEST_H
              #define CTEST_H
              #include '
  siapoplib <- paste('"', .libPaths(), "/siapopr/include/constantGlobalStructs.h", '"', sep = "")
  header2 <- '
    #include <gsl/gsl_randist.h>
    #ifdef __cplusplus
    extern "C" {
      #endif
      void customdist(double* fitness, struct FitnessParameters *fit_params, gsl_rng* rng);
      #ifdef __cplusplus
    }
    #endif
    #endif'
  writeLines(paste(header1, siapoplib, header2, sep = ""), con, sep = "")
  close(con)

  con <- file(cppfile)
  headerfile <- paste(cpproot, ".h", sep = "")
  include_header <- paste('#include "', .pop(headerfile, "/"), '"', sep = "")
  cppdat <- readLines(con)
  if(!(include_header %in% cppdat)) writeLines(c(include_header, cppdat), con, sep = "\n")
  close(con)

  concopy <- file(paste(cppfile, ".backup", sep = ""))
  writeLines(cppdat, concopy)
  close(concopy)

  compile <- paste("R CMD COMPILE ", cpproot, ".", cppsuffix, sep = "")
  # Need to make windows version
  shlib <- paste("R CMD SHLIB -o ", cpproot, ".so ",  cpproot, ".o -dynamiclib -lgsl", sep = "")
  system(compile)
  system(shlib)
}



##------------------------------------------------------------------------
#' compile_custom_newclone
#'
#' Allows the user to create a custom plugin in C++ for use in
#' SIApop. The plugin would allow new inheritance scenarios for a new clone
#' enterring the population.
#'
#' @param cppfile a character string giving the path name of a cpp file
#'
#' @export
#' @examples
#' \dontrun{
#' compile_custom_newclone(cppfile = "./plugin.cpp")
#' }
compile_custom_newclone <- function(cppfile){
  cpproot <- .pop_off(cppfile, ".", fixed = T)
  cppsuffix <- .pop(cppfile, ".", fixed = T)
  con <- file(paste(cpproot, ".h", sep = ""))

  header1 <- '#ifndef CCLONE_H
#define CCLONE_H'
  structurelib <- paste('#include "', .libPaths(), "/siapopr/include/constantGlobalStructs.h", '"', sep = "")
  functionlib <- paste('#include "', .libPaths(), "/siapopr/include/constantRVFunctions.h", '"', sep = "")
  header2 <- '
#include <gsl/gsl_randist.h>
#include <math.h>

  #ifdef __cplusplus
  extern "C" {
  #endif
  void customclone(struct clone *new_clone, struct clone *parent_clone,
  struct FitnessParameters* fit_params, struct MutationParameters* mut_params,
  struct PunctuationParameters* punct_params,
  struct EpistaticParameters* epi_params, gsl_rng* rng,
  void (*ConstantGenerateFitness)(double *, struct FitnessParameters*, gsl_rng*));
  #ifdef __cplusplus
  }
  #endif
  #endif'
  writeLines(paste(header1, structurelib, functionlib, header2, sep = "\n"), con, sep = "")
  close(con)

  con <- file(cppfile)
  headerfile <- paste(cpproot, ".h", sep = "")
  include_header <- paste('#include "', .pop(headerfile, "/"), '"', sep = "")
  cppdat <- readLines(con)
  if(!(include_header %in% cppdat)) writeLines(c(include_header, cppdat), con, sep = "\n")
  close(con)

  concopy <- file(paste(cppfile, ".backup", sep = ""))
  writeLines(cppdat, concopy)
  close(concopy)

  compile <- paste("R CMD COMPILE ", cpproot, ".", cppsuffix, sep = "")
  # Need to make windows version
  shlib <- paste("R CMD SHLIB -o ", cpproot, ".so ",  cpproot, ".o -dynamiclib -lgsl", sep = "")
  system(compile)
  system(shlib)
}

