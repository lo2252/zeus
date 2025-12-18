#' @title Importing raw erg.abf files
#'
#' @description
#' This function will import the raw .abf and convert it into a usable data.frame for manipulation.\n
#' This utilizes the package 'readABF', please see CRAN publications for this package for more information.
#'
#' @author Logan Ouellette
#'
#'
#'
#'

import_erg <- function(erg_file_name) {

erg_file <- readABF::as.data.frame.ABF()

}











