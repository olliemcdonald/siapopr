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
  if(cppfile[length(cppfile)] != "cpp")
  {
    hfile <- c(cppfile, "h")
    cppfile <- c(cppfile, "cpp")
  }
  else
  {
    hfile <- c(cppfile[1:(length(cppfile)-1)], "h")
  }

  cppfile <- paste(cppfile, collapse = ".")
  hfile <- paste(hfile, collapse = ".")


  cpp_location <- paste(.libPaths()[1], "/siapopr/extras/customdist_template.cpp", sep = "")
  h_location <- paste(.libPaths()[1], "/siapopr/extras/customdist_template.h", sep = "")

  file.copy(cpp_location, cppfile)
  file.copy(h_location, hfile)

  # Update cpp file's include statement
  cpp_con <- file(cppfile)
  cpp_lines <- readLines(cpp_con)
  include_header <- paste("#include \"", hfile, "\"", sep = "")
  cpp_lines <- c(include_header, cpp_lines)
  write(cpp_lines, cppfile)
  close(cpp_con)
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
  if(cppfile[length(cppfile)] != "cpp")
  {
    hfile <- c(cppfile, "h")
    cppfile <- c(cppfile, "cpp")
  }
  else
  {
    hfile <- c(cppfile[1:(length(cppfile)-1)], "h")
  }
  cppfile <- paste(cppfile, collapse = ".")
  hfile <- paste(hfile, collapse = ".")
  cpp_location <- paste(.libPaths()[1], "/siapopr/extras/custom_newclone_template.cpp", sep = "")
  h_location <- paste(.libPaths()[1], "/siapopr/extras/custom_newclone_template.h", sep = "")

  file.copy(cpp_location, cppfile)
  file.copy(h_location, hfile)

  # Update header file's include statements
  h_con <- file(hfile)
  h_lines <- readLines(h_con)
  h_lines[3] <- paste('#include "', .libPaths()[1], "/siapopr/include/constantGlobalStructs.h", '"', sep = "")
  h_lines[4] <- paste('#include "', .libPaths()[1], "/siapopr/include/constantRVFunctions.h", '"', sep = "")
  write(h_lines, hfile)
  close(h_con)

  # Update cpp file's include statement
  cpp_con <- file(cppfile)
  cpp_lines <- readLines(cpp_con)
  include_header <- paste("#include \"", hfile, "\"", sep = "")
  cpp_lines <- c(include_header, cpp_lines)
  write(cpp_lines, cppfile)
  close(cpp_con)
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

  concopy <- file(paste(cppfile, ".backup", sep = ""))
  writeLines(cppfile, concopy)
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

  concopy <- paste(cppfile, ".backup", sep = "")
  file.copy(cppfile, concopy)
  close(concopy)

  compile <- paste("R CMD COMPILE ", cpproot, ".", cppsuffix, sep = "")
  # Need to make windows version
  shlib <- paste("R CMD SHLIB -o ", cpproot, ".so ",  cpproot, ".o -dynamiclib -lgsl", sep = "")
  system(compile)
  system(shlib)
}

