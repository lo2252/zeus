#' @title Importing raw erg.abf files
#'
#' @description
#' This function will import the raw .abf and convert it into a usable data.frame for manipulation.\n
#' This utilizes the package 'readABF', please see CRAN publications for this package for more information.
#'
#' @author Logan Ouellette


#' @title Read an ABF file and tryutn a tidy data frame
#' @param path Path to a .abf file
#' @param ... Passed to readABF::readABF()
#' @return A data.frame
#' @export

read_abf_to_df <- function(path, ...) {

  if(!file.exists(path)) stop("File not found: ", path, call. = FALSE)

  if (!requireNamespace("readABF", quietly = TRUE)) {
    stop("Package 'readABF' is required. Install it with install.packages('readABF').", call. = FALSE)
  }

  # Reads raw ABF file
  raw_abf <- readABF::readABF(path, ...)
  # Calls df conversion function
  abf_as_df(raw_abf)


}

#' @title Convert readABF output to a standardizes data frame
#'
#'
#'
#'

abf_as_df <- function(raw_abf) {
  df <- as.data.frame(raw_abf, sweep = 1)

}










