#' Importing raw erg.abf files
#'
#' @description
#' This function imports a raw .abf file and stores it as a `zeus_abf_raw`
#' object for downstream conversion and protocol annotation.
#'
#' This utilizes the package `readABF`; please see CRAN documentation for
#' more information.
#'
#' @author Logan Ouellette


#' @title Read an ABF file into zeus raw object
#' @param path Path to a .abf file
#' @param ... Passed to readABF::readABF()
#' @return An object of class `zeus_abf_raw` wrapping the readABF output
#' @export
read_abf_raw <- function(path, ...) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)

  if (!requireNamespace("readABF", quietly = TRUE)) {
    stop(
      "Package 'readABF' is required. Install it with install.packages('readABF').",
      call. = FALSE
    )
  }

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
#' @param protocol Character string specifying the protocol. Supported:
#'   `"default"`, `"C1"`, `"C0"`.
#' @param nd_start,nd_end,nd_step,repeats_per_level,n_protocol_repeats Protocol
#'   settings used by `"default"` and `"C1"`. `"C0"` uses a fixed custom sequence.
#'
#' @return A tibble/data.frame in standardized long format.
#' @export
read_abf <- function(
    x,
    ...,
    add_stim = TRUE,
    calib = NULL,
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  protocol <- match.arg(protocol)

  raw_abf <- NULL

  if (inherits(x, "zeus_abf_raw")) {
    raw_abf <- x$raw
  } else if (is.character(x) && length(x) == 1) {
    raw_abf <- read_abf_raw(x, ...)$raw
  } else {
    stop("'x' must be a file path or a 'zeus_abf_raw' object.", call. = FALSE)
  }

  df_long <- abf_as_df_long(raw_abf)

  if (isTRUE(add_stim)) {
    df_long <- add_stimulus_cols_protocol(
      df_long = df_long,
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

# Allowed Treatment values -----------------------------------------------------
#' Allowed treatment group values
#'
#' Internal helper returning the controlled vocabulary of supported
#' treatment group identifiers used in ZEUS metadata.
#'
#' @return A character vector of allowed treatment group labels.
#' @keywords internal


.allowed_treatment_groups <- function() {
  c("SYS Water",
    "BPA",
    "DMSO",
    "1-850",
    "T3",
    "ICI",
    "EE2",
    "BPA/ICI",
    "BPA/1-850",
    "BPA/1-850/ICI",
    "user_input")
}

#' Allowed ERG age categories
#'
#' Internal helper returning the allowed age categories for ERG recordings.
#'
#' @return A character vector of allowed ERG age values.
#' @keywords internal


.allowed_erg_age <- function() {
  c("Larval", "Adult")
}

# Convert to long DF -----------------------------------------------------------

#' Convert readABF output to a standardized long data frame
#'
#' @param raw_abf An object of class `ABF`.
#' @return A tibble in long format with sweep, time, channel, units, and value.
#' @keywords internal
abf_as_df_long <- function(raw_abf) {
  stopifnot(inherits(raw_abf, "ABF"))

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
    time <- (seq_len(npts) - 1) * dt

    tibble::tibble(
      sweep = s,
      time = rep(time, times = ncol(m)),
      channel = rep(ch_names, each = npts),
      units = rep(ch_units, each = npts),
      value = as.vector(m)
    )
  })
}


# Convert to wide DF -----------------------------------------------------------

#' Convert readABF output to a standardized wide data frame
#'
#' @param raw_abf An object of class `ABF`.
#' @return A tibble in wide format with one row per time point per sweep.
#' @keywords internal
abf_as_df_wide <- function(raw_abf) {
  stopifnot(inherits(raw_abf, "ABF"))

  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required for this method.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for this method.", call. = FALSE)
  }

  dt <- raw_abf$samplingIntervalInSec
  ch_names <- raw_abf$channelNames

  purrr::imap_dfr(raw_abf$data, function(m, s) {
    npts <- nrow(m)
    time <- (seq_len(npts) - 1) * dt

    df <- tibble::as_tibble(m, .name_repair = "minimal")
    names(df) <- ch_names

    tibble::add_column(
      df,
      sweep = s,
      time = time,
      .before = 1
    )
  })
}


# Calibration table ------------------------------------------------------------

#' Default Neutral Density (ND) to irradiance calibration table
#'
#' @keywords internal
.default_stim_calib <- function() {
  data.frame(
    stim_nd = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0),
    stim_irradiance_log10 = c(
      -5.977, -5.321, -5.003,
      -4.491, -3.957, -3.457, -2.927
    )
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


# Protocol builders ------------------------------------------------------------

#' Build the default ND-only protocol
#'
#' @param nd_start Starting ND level.
#' @param nd_end Ending ND level.
#' @param nd_step Step size between ND levels.
#' @param repeats_per_level Number of repeated sweeps per ND level.
#' @param n_protocol_repeats Number of full protocol repeats.
#'
#' @return A data.frame describing the sweep protocol.
#' @keywords internal
make_protocol_default <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  nd_levels <- seq(nd_start, nd_end, by = nd_step)
  one_protocol <- rep(nd_levels, each = repeats_per_level)
  stim_nd <- rep(one_protocol, times = n_protocol_repeats)

  data.frame(
    stim_order = seq_along(stim_nd),
    stim_type = rep("flash", length(stim_nd)),
    stim_nd = stim_nd,
    wavelength = NA_real_,
    stringsAsFactors = FALSE
  )
}


#' Build the C1 white-light protocol
#'
#' @param nd_start Starting ND level.
#' @param nd_end Ending ND level.
#' @param nd_step Step size between ND levels.
#' @param repeats_per_level Number of repeated sweeps per ND level.
#' @param n_protocol_repeats Number of full protocol repeats.
#'
#' @return A data.frame describing the sweep protocol.
#' @keywords internal
make_protocol_C1 <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  nd_levels <- seq(nd_start, nd_end, by = nd_step)
  one_protocol <- rep(nd_levels, each = repeats_per_level)
  stim_nd <- rep(one_protocol, times = n_protocol_repeats)

  data.frame(
    stim_order = seq_along(stim_nd),
    stim_type = rep("flash", length(stim_nd)),
    stim_nd = stim_nd,
    wavelength = rep(570, length(stim_nd)),
    stringsAsFactors = FALSE
  )
}


#' Build the C0 spectral protocol
#'
#' In the spectral protocol, each wavelength has its own ordered ND sequence.
#' Each ND value is repeated 4 times before advancing to the next wavelength.
#' The full spectral schedule contains 10 wavelength blocks x 7 ND values x
#' 4 repeats = 280 sweeps.
#'
#' @param repeats_per_level Number of repeats for each ND value within a wavelength.
#'
#' @return A data.frame describing the sweep protocol.
#' @keywords internal
make_protocol_C0 <- function(repeats_per_level = 4) {
  spectral_blocks <- list(
    list(wavelength = 650, stim_nd = c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0)),
    list(wavelength = 570, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 490, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 410, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 330, stim_nd = c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0)),
    list(wavelength = 650, stim_nd = c(4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5)),
    list(wavelength = 610, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 530, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 450, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    list(wavelength = 370, stim_nd = c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0))
  )

  protocol_df <- do.call(
    rbind,
    lapply(spectral_blocks, function(block) {
      nd_seq <- rep(block$stim_nd, each = repeats_per_level)

      data.frame(
        stim_type = rep("flash", length(nd_seq)),
        stim_nd = nd_seq,
        wavelength = rep(block$wavelength, length(nd_seq)),
        stringsAsFactors = FALSE
      )
    })
  )

  protocol_df$stim_order <- seq_len(nrow(protocol_df))
  protocol_df <- protocol_df[, c("stim_order", "stim_type", "stim_nd", "wavelength")]

  rownames(protocol_df) <- NULL
  protocol_df
}


#' Build protocol table based on selected protocol
#'
#' @param protocol Character string specifying protocol.
#' @param nd_start Starting ND level.
#' @param nd_end Ending ND level.
#' @param nd_step Step size between ND levels.
#' @param repeats_per_level Number of repeated sweeps per ND level.
#' @param n_protocol_repeats Number of full protocol repeats.
#'
#' @return A protocol table with `stim_order`, `stim_type`, `stim_nd`, and `wavelength`.
#' @keywords internal
make_protocol_table <- function(
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10
) {
  protocol <- match.arg(protocol)

  switch(
    protocol,
    default = make_protocol_default(
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats
    ),
    C1 = make_protocol_C1(
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats
    ),
    C0 = make_protocol_C0(
      repeats_per_level = repeats_per_level
    )
  )
}


# Protocol joiner --------------------------------------------------------------

#' Add protocol-derived stimulus columns to long-form ABF data
#'
#' @param df_long Long-form ABF data.
#' @param calib Optional calibration table.
#' @param protocol Character string specifying protocol.
#' @param nd_start Starting ND level for default/C1 protocols.
#' @param nd_end Ending ND level for default/C1 protocols.
#' @param nd_step Step size for default/C1 protocols.
#' @param repeats_per_level Number of repeated sweeps per ND level.
#' @param n_protocol_repeats Number of full protocol repeats.
#'
#' @return The input data frame with stimulus metadata joined by sweep.
#' @keywords internal
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
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

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
    stim_order = protocol_tbl$stim_order[seq_len(n_map)],
    stim_type = protocol_tbl$stim_type[seq_len(n_map)],
    stim_nd = protocol_tbl$stim_nd[seq_len(n_map)],
    wavelength = protocol_tbl$wavelength[seq_len(n_map)]
  ) |>
    dplyr::mutate(
      stim_irradiance_log10 = nd_to_irradiance_log10(stim_nd, calib = calib)
    )

  dplyr::left_join(df_long, stim_tbl, by = "sweep")
}
