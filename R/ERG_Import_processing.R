#' Importing raw erg.abf files
#'
#' @description
#' This function will import the raw.abf and convert it into a usable data.frame for manipulation.\n
#' This utilizes the package 'readABF', please see CRAN publications for this package for more information.
#'
#' @author Logan Ouellette


#' @title Read an ABF file into zeus raw object
#' @param path Path to a .abf file
#' @param ... Passed to readABF::readABF()
#' @return An object of class 'zeus_abf_raw' wrapping the readABF output
#' @export
read_abf_raw <- function(path, ...) {
  if(!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (!requireNamespace("readABF", quietly = TRUE)) {
    stop("Package 'readABF' is required. Install it with install.packages('readABF').", call. = FALSE)
  }

  # Reads raw ABF file
  raw_abf <- readABF::readABF(path, ...)

  structure(
    list(
      raw = raw_abf,
      path = normalizePath(path, winslash = "/", mustWork = FALSE)
    ),
    class = "zeus_abf_raw"
  )
}

# Converting imported raw file into usable data.frame --------------------------

#' Read an ABF file and return a standardized long data frame
#'
#' @param x Either a file path to a .abf, or a `zeus_abf_raw` object.
#' @param ... Passed to readABF::readABF() only when `x` is a file path.
#' @param add_stim Logical; attach protocol stimulus columns?
#' @param calib Optional calibration table passed to nd_to_irradiance_log10().
#' @param nd_start,nd_end,nd_step,repeats_per_level,n_protocol_repeats Protocol settings.
#'
#' @return A tibble/data.frame in standardized long format.
#' @export

read_abf <- function(
  x,
  ...,
  add_stim = TRUE,
  calib = NULL,
  protocol = (c("default", "C1", "C0")),
  nd_start = 3.0,
  nd_end = 6.0,
  nd_step = 0.5,
  repeats_per_level = 4,
  n_protocol_repeats = 10) {

  # Listing Protocol
  protocol <- match.arg(protocol)

  # Accept either a path or a zeus_abf_raw object
  raw_abf <- NULL

  if (inherits(x, "zeus_abf_raw")) {
    raw_abf <- x$raw
  } else if (is.character(x) && length(x) == 1) {
    raw_abf <- read_abf_raw(x, ...)$raw
  } else {
    stop("'x' must be a file path or a 'zeus_abf_raw' object.", call = FALSE)
  }

  # Convert to long
  df_long <- abf_as_df_long(raw_abf)

  # Optional adding stimulus, default TRUE
  if (isTRUE(add_stim)) {
    df_long <- add_stimulus_cols_protocol(
      df_long,
      calib = calib,
      protocol = protocol,
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats
    )
  }
  df_long
}



# Convert to long DF -----------------------------------------------------------

#' Convert readABF output to a standardizes data frame
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

# Default Protocol:
# ND 3.0 to 6.0 by 0.5, each ND repeats 4 times, repeated 10 times
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

# C1 (White light):
# same ND protocol as above, but wavelength is fixed at 570
make_protocol_C1 <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  stim_nd <- make_stim_nd_sweep_protocol(
    nd_start = nd_start,
    nd_end = nd_end,
    nd_step = nd_step,
    repeats_per_level = repeats_per_level,
    n_protocol_repeats = n_protocol_repeats
  )

  data.frame(
    stim_nd = stim_nd,
    wavelength = 570
  )
}


# C0 (Spectral):
# wavelength and ND
make_protocol_C0 <- function() {

  spectral_blocks <- list(
    list(wavelength = 650, stim_nd = c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, NA)),
    list(wavelength = 570, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 490, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 410, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 330, stim_nd = c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, NA)),
    list(wavelength = 650, stim_nd = c(4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, NA)),
    list(wavelength = 610, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 530, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 450, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA)),
    list(wavelength = 370, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, NA))
  )

  out <- do.call(
    rbind,
    lapply(spectral_blocks, function(block) {
      data.frame(
        stim_nd = block$stim_nd,
        wavelength = rep(block$wavelength, length(block$stim_nd))
      )
    })
  )

  rownames(out) <- NULL
  out
}


# Helper to build protocol table
make_protocol_table <- function(
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  protocol <- match.arg(protocol)

  if (protocol == "default") {
    data.frame(
      stim_nd = make_stim_nd_sweep_protocol(
        nd_start = nd_start,
        nd_end = nd_end,
        nd_step = nd_step,
        repeats_per_level = repeats_per_level,
        n_protocol_repeats = n_protocol_repeats
      ),
      wavelength = NA_real_
    )
  } else if (protocol == "C1") {
    make_protocol_C1(
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats
    )
  } else if (protocol == "C0") {
    make_protocol_C0()
  } else {
    stop("Unsupported protocol: ", protocol, call. = FALSE)
  }
}

# Adds stims to protocol
add_stimulus_cols_protocol <- function(
    df_long,
    calib = NULL,
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.", call. = FALSE)
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Need tibble.", call. = FALSE)

  protocol <- match.arg(protocol)

  protocol_tbl <- make_protocol_table(
    protocol = protocol,
    nd_start = nd_start,
    nd_end = nd_end,
    nd_step = nd_step,
    repeats_per_level = repeats_per_level,
    n_protocol_repeats = n_protocol_repeats
  )

  sweeps_present <- sort(unique(df_long$sweep))
  n_map <- min(length(sweeps_present), nrow(protocol_tbl))

  stim_tbl <- tibble::tibble(
    sweep = sweeps_present[seq_len(n_map)],
    stim_nd = protocol_tbl$stim_nd[seq_len(n_map)],
    wavelength = protocol_tbl$wavelength[seq_len(n_map)]
  ) |>
    dplyr::mutate(
      stim_irradiance_log10 = dplyr::if_else(
        is.na(stim_nd),
        NA_real_,
        nd_to_irradiance_log10(stim_nd, calib = calib)
      )
    )

  dplyr::left_join(df_long, stim_tbl, by = "sweep")
}


