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
