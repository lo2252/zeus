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
  abf_as_df_wide(raw_abf) # Change to either _long or _wide


}

#' @title Convert readABF output to a standardizes data frame
#'
#'
#'
#'

# Function in long format
abf_as_df_long <- function(raw_abf) {
  stopifnot(inherits(raw_abf, "ABF"))

  # Stops if missing packages
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required for this method.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for this method.", call. = FALSE)
  }

  dt <- raw_abf$samplingIntervalInSec
  ch_names <- raw_abf$channelNames
  ch_units <- raw_abf$channelUnits

  purrr::imap_dfr(raw_abf$data, function(m, s) {
    npts <- nrow(m)
    time <- (seq_len(npts) - 1) * dt # seconds

    # Creating long tibble
    tibble::tibble(
      sweep = s,
      time = rep(time, times = ncol(m)),
      channel = rep(ch_names, each = npts),
      units = rep(ch_units, each = npts),
      value = as.vector(m)
    )
  })

  }


# Function in wide
abf_as_df_wide <- function(raw_abf) {
  stopifnot(inherits(raw_abf, "ABF"))

  dt <- raw_abf$samplingIntervalInSec
  ch_names <- raw_abf$channelNames

  purrr::imap_dfr(raw_abf$data, function(m,s) {
    npts <- nrow(m)
    time <- (seq_len(npts) - 1)* dt

    df <- tibble::as_tibble(m, .name_repair = "minimal")
    names(df) <- ch_names
    tibble::add_column(df,
                       sweep = s,
                       time = time,
                       .before = 1)
  })
}





