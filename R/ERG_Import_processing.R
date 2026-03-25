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
#' @param treatment_group Treatment label. One of `"SYS Water"`, `"BPA"`,
#'   `"DMSO"`, `"1-850"`, `"T3"`, `"ICI"`, `"EE2"`, `"BPA/ICI"`,
#'   `"BPA/1-850"`, `"BPA/1-850/ICI"`, or `"user_input"`.
#' @param treatment_group_custom Character string used when
#'   `treatment_group = "user_input"`.
#' @param date_of_fertilization Date of fertilization for the fish.
#' @param erg_age Age at ERG, one of `"Larval"` or `"Adult"`.
#' @param apply_boxcar Logical; if `TRUE`, applies a boxcar filter to selected
#'   response channels before stimulus annotation.
#' @param boxcar_channel_pattern Character string used to identify the channel
#'   to smooth. Default is `"DAM80"`.
#' @param boxcar_k Optional integer window size for the boxcar filter. If not
#'   supplied, `boxcar_width` is used to determine the window size.
#' @param boxcar_width Optional numeric smoothing width in the same units as
#'   the `time` column. Used only when `boxcar_k` is `NULL`.
#' @param keep_raw_boxcar Logical; if `TRUE`, stores the original unsmoothed
#'   signal in a new column called `value_raw`.
#'
#' @return A tibble/data.frame in standardized long format.
#' @export
zeus_import <- function(
    x,
    ...,
    add_stim = TRUE,
    calib = NULL,
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10,
    treatment_group = NULL,
    treatment_group_custom = NULL,
    date_of_fertilization = NA,
    erg_age = NULL,
    apply_boxcar = FALSE,
    boxcar_channel_pattern = "DAM80",
    boxcar_k = NULL,
    boxcar_width = NULL,
    keep_raw_boxcar = FALSE
) {

  protocol <- match.arg(protocol)

  if (inherits(x, "zeus_abf_raw")) {
    raw_abf <- x$raw
  } else if (is.character(x) && length(x) == 1) {
    raw_abf <- read_abf_raw(x, ...)$raw
  } else {
    stop("'x' must be a file path or a 'zeus_abf_raw' object.", call. = FALSE)
  }

  df_long <- abf_as_df_long(raw_abf)

  if (isTRUE(apply_boxcar)) {

    if (is.null(boxcar_k) && is.null(boxcar_width)) {
      boxcar_width <- 0.002
    }

    df_long <- zeus_boxcar_filter(
      df_long = df_long,
      channel_pattern = boxcar_channel_pattern,
      k = boxcar_k,
      boxcar_width = boxcar_width,
      keep_raw = keep_raw_boxcar
    )
  }

  if (isTRUE(add_stim)) {
    df_long <- add_stimulus_cols_protocol(
      df_long = df_long,
      calib = calib,
      protocol = protocol,
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats,
      treatment_group = treatment_group,
      treatment_group_custom = treatment_group_custom,
      date_of_fertilization = date_of_fertilization,
      erg_age = erg_age
    )
  }

  df_long
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
    stop("Package 'purrr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
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


# Boxcar filter -----------------------------------------------------------

#' Apply a boxcar filter to selected channels in long-format ABF data
#'
#' Applies a centered moving-average (boxcar) filter to the `value` column
#' within each sweep for channels matching `channel_pattern`.
#'
#' @param df_long Long-format ABF data.
#' @param channel_pattern Character string identifying channel to smooth.
#' @param k Optional integer window size for the filter. If `NULL`, `k` is
#'   computed from `boxcar_width` and the spacing in `time`.
#' @param boxcar_width Optional numeric smoothing width in the same units as
#'   `df_long$time`. Used only when `k` is `NULL`.
#' @param keep_raw Logical; if `TRUE`, store original signal in `value_raw`.
#'
#' @return Data frame with smoothed `value`.
#' @keywords internal
zeus_boxcar_filter <- function(df_long,
                               channel_pattern = "DAM80",
                               k = NULL,
                               boxcar_width = NULL,
                               keep_raw = FALSE) {

  required_cols <- c("sweep", "time", "channel", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(k)) {
    if (is.null(boxcar_width) || !is.numeric(boxcar_width) ||
        length(boxcar_width) != 1 || is.na(boxcar_width) || boxcar_width <= 0) {
      stop(
        "When `k` is NULL, `boxcar_width` must be a single numeric value > 0.",
        call. = FALSE
      )
    }

    unique_time <- sort(unique(df_long$time))

    if (length(unique_time) < 2) {
      stop(
        "Not enough unique time points to estimate sampling interval.",
        call. = FALSE
      )
    }

    dt <- stats::median(diff(unique_time), na.rm = TRUE)

    if (!is.finite(dt) || dt <= 0) {
      stop(
        "Could not determine a valid sampling interval from `time`.",
        call. = FALSE
      )
    }

    k <- as.integer(round(boxcar_width / dt))
    k <- max(1L, k)

    if (k %% 2 == 0) {
      k <- k + 1L
    }
  } else {
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1) {
      stop("`k` must be a single integer >= 1.", call. = FALSE)
    }

    k <- as.integer(k)

    if (k %% 2 == 0) {
      k <- k + 1L
    }
  }

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_long)) {
    df_long$value_raw <- df_long$value
  }

  idx <- grepl(channel_pattern, df_long$channel)

  if (!any(idx)) {
    return(df_long)
  }

  df_target <- df_long[idx, , drop = FALSE]
  df_other  <- df_long[!idx, , drop = FALSE]

  split_groups <- split(
    df_target,
    interaction(df_target$sweep, df_target$channel, drop = TRUE)
  )

  filtered_groups <- lapply(split_groups, function(dat) {
    dat <- dat[order(dat$time), , drop = FALSE]

    filt <- stats::filter(
      dat$value,
      rep(1 / k, k),
      sides = 2
    )

    filt <- as.numeric(filt)

    na_idx <- is.na(filt)
    filt[na_idx] <- dat$value[na_idx]

    dat$value <- filt
    dat
  })

  out <- do.call(rbind, c(filtered_groups, list(df_other)))
  out <- out[order(out$sweep, out$channel, out$time), , drop = FALSE]
  rownames(out) <- NULL

  out
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


# Treatment and Age Metadata Helper --------------------------------------------
#' Construct sample-level metadata for ERG recordings
#'
#' Internal helper used to validate and construct meta data describing
#' the sample associated with the ERG recording. This includes the treatment
#' group identity, date of fertilization, and developmental stage at the time
#' of ERG recording.
#'
#' Treatment groups are confirmed against a controlled list by
#' `.allowed_treatment_grous()`. If `"user_input"` is specified, a custom
#' treatment label must be provided via `treatment_group_custom`.
#'
#' ERG age is validated against `.allowed_erg_age()`.
#'
#' @param treatment_group Character string specifying treatment group.
#'  Must be one of: `"SYS Water"`, `"BPA"`, `"DMSO"`, `"1-850"`, `"T3"`,
#'   `"ICI"`, `"EE2"`, `"BPA/ICI"`, `"BPA/1-850"`, `"BPA/1-850/ICI"`,
#'   or `"user_input"`.
#'
#' @param treatment_group_custom Character string used when `treatment_group`
#'   = `"user_input"` to specify a custom treatment label.
#'
#' @param date_of_fertilization Date of fertilization for the fish. Accepts
#' `Date` objects or date-like value, using `as.Date()` to coerce
#' for the format.
#'
#' @param erg_date Character string describing developmental stage at the time
#' of ERG recording. Must be `"larval"` or `"adult"`.
#'
#' @return A data frame containing validated sample metadata with columns:
#' \describe{
#'  \item{treatment_group}{Treatment group label.}
#'  \item{date_of_fertilization}{Date of fertilization for the sample.}
#'  \item{erg_age}{Developmental stage at ERG recording.}
#'  }
#' @importFrom lubridate mdy
#' @keywords internal

make_sample_metadata <- function(
    treatment_group = NULL,
    treatment_group_custom = NULL,
    date_of_fertilization = NA,
    erg_age = NULL
) {
  allowed_treatments <- .allowed_treatment_groups()
  allowed_ages <- .allowed_erg_age()

  if (is.null(treatment_group)) {
    treatment_group <- NA_character_
  }

  if (!is.na(treatment_group) && !treatment_group %in% allowed_treatments) {
    stop(
      "'treatment_group' must be one of: ",
      paste(allowed_treatments, collapse = ", "),
      call. = FALSE
    )
  }

  if (identical(treatment_group, "user_input")) {
    if (is.null(treatment_group_custom) || !nzchar(trimws(treatment_group_custom))) {
      stop(
        "When treatment_group = 'user_input', you must supply ",
        "'treatment_group_custom'.",
        call. = FALSE
      )
    }
    treatment_group_final <- trimws(treatment_group_custom)
  } else {
    treatment_group_final <- treatment_group
  }

  if (is.null(erg_age)) {
    erg_age <- NA_character_
  }

  if (!is.na(erg_age) && !erg_age %in% allowed_ages) {
    stop(
      "'erg_age' must be one of: ",
      paste(allowed_ages, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(date_of_fertilization) || length(date_of_fertilization) == 0) {
    date_of_fertilization <- as.Date(NA)
  } else if (inherits(date_of_fertilization, "Date")) {
    # leave as-is
  } else if (is.character(date_of_fertilization)) {
    date_of_fertilization <- lubridate::mdy(date_of_fertilization, quiet = TRUE)

    if (is.na(date_of_fertilization)) {
      stop(
        "'date_of_fertilization' must be a valid date in month/day/year format ",
        "(e.g., '03/23/2026') or a Date object.",
        call. = FALSE
      )
    }
  } else if (is.na(date_of_fertilization)) {
    date_of_fertilization <- as.Date(NA)
  } else {
    stop(
      "'date_of_fertilization' must be NA, a Date object, or a character string ",
      "in month/day/year format.",
      call. = FALSE
    )
  }

  data.frame(
    treatment_group = treatment_group_final,
    date_of_fertilization = date_of_fertilization,
    erg_age = erg_age,
    stringsAsFactors = FALSE
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
    n_protocol_repeats = 10,
    treatment_group = NULL,
    treatment_group_custom = NULL,
    date_of_fertilization = NA,
    erg_age = NULL
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
    ) |>
    dplyr::distinct(sweep, .keep_all = TRUE)

  if (any(is.na(df_long$sweep))) {
    stop("`df_long$sweep` contains NA values; cannot join safely.", call. = FALSE)
  }

  if (any(is.na(stim_tbl$sweep))) {
    stop("`stim_tbl$sweep` contains NA values; cannot join safely.", call. = FALSE)
  }

  if (anyDuplicated(stim_tbl$sweep)) {
    dupes <- stim_tbl$sweep[duplicated(stim_tbl$sweep)]
    stop(
      "Duplicate sweep values found in `stim_tbl`: ",
      paste(unique(dupes), collapse = ", "),
      call. = FALSE
    )
  }

  sample_meta <- make_sample_metadata(
    treatment_group = treatment_group,
    treatment_group_custom = treatment_group_custom,
    date_of_fertilization = date_of_fertilization,
    erg_age = erg_age
  )

  df_long |>
    dplyr::mutate(sweep = as.character(sweep)) |>
    dplyr::left_join(
      stim_tbl |> dplyr::mutate(sweep = as.character(sweep)),
      by = "sweep",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      treatment_group = sample_meta$treatment_group,
      date_of_fertilization = sample_meta$date_of_fertilization,
      erg_age = sample_meta$erg_age
    )
}
