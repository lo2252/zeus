#' @title Importing raw erg.abf files
#'
#' @description
#' This function will import the raw .abf and convert it into a usable data.frame for manipulation.\n
#' This utilizes the package 'readABF', please see CRAN publications for this package for more information.
#'
#' @author Logan Ouellette


#' @title Read an ABF file and convert to a tidy data frame
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
  abf_as_df_long(raw_abf) # Change to either _long or _wide


}

# Convert to long DF -----------------------------------------------------------

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



# Convert to wide DF ------------------------------------------------------


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


# Adding stimulus to column -----------------------------------------------

#' Default Neutral Density (ND), irradiance calibration table
#'
#' @keywords internal

.default_stim_calib <- function() {
  data.frame(
    stim_nd = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0),
    stim_irradiance_log10 = c(-5.977, -5.321, -5.003,
                               -4.491, -3.957, -3.457, -2.927)
  )
}

#' Convert ND to log10 irradiance
#'
#' @param stim_nd Numeric vector of ND values.
#' @param calib Optional calibration data frame with columns
#'   `stim_nd` and `stim_irradiance_log10`.
#'
#' @return Numeric vector of log10 irradiance values.
#' @export

nd_to_irradiance_log10 <- function(stim_nd, calib = NULL) {
  if (is.null(calib)) calib <- .default_stim_calib()

  x <- as.numeric(unlist(calib[["stim_nd"]]))
  y <- as.numeric(unlist(calib[["stim_irradiance_log10"]]))

  o <- order(x)
  stats::approx(x = x[o], y = y[o], xout = stim_nd, rule = 2)$y
}


# Stimulus Schedule and irradiance mapping --------------------------------

# Protocol: ND 3.0 to 6.0 by 0.5, each ND repeats 4 times, repeating the full sequence
# 10 time = 280 sweeps
make_stim_nd_sweep_protocol <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  nd_levels <- seq(nd_start, nd_end, by = nd_step)
  one_protocol <- rep(nd_levels, each = repeats_per_level)
  full <- rep(one_protocol, times = n_protocol_repeats)
  full
}

add_stimulus_cols_protocol <- function(
    df_long,
    calib = NULL,
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.", call. = FALSE)
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Need tibble.", call. = FALSE)

  # Building protocol for ND schedule (length 280)
  stim_nd_seq <- make_stim_nd_sweep_protocol(
    nd_start = nd_start,
    nd_end = nd_end,
    nd_step = nd_step,
    repeats_per_level = repeats_per_level,
    n_protocol_repeats = n_protocol_repeats
  )

  # Map protocol onto sweeps existing in data (keeping extra sweeps)
  sweeps_present <- sort(unique(df_long$sweep))
  n_map <- min(length(sweeps_present), length(stim_nd_seq))

  stim_tbl <- tibble::tibble(
    sweep = sweeps_present[seq_len(n_map)],
    stim_nd = stim_nd_seq[seq_len(n_map)]) |>

    dplyr::mutate(
      stim_irradiance_log10 = nd_to_irradiance_log10(stim_nd, calib = calib)
    )

  # Left join keeps all rows, any sweeps outside of protocol recieve NA

  dplyr::left_join(df_long, stim_tbl, by = "sweep")

}



