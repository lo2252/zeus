# Importing Raw ERG -------------------------------------------------------

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

# Import/Clean Data -------------------------------------------------------

#' @description
#' Imports Axon Binary Format (ABF) electrophysiology data and converts it into
#' a standardized long-format data frame suitable for downstream analysis within
#' ZEUS.
#'
#' The import workflow can:
#' \enumerate{
#'   \item read raw ABF data,
#'   \item convert the file to ZEUS long format,
#'   \item append stimulus metadata,
#'   \item assign protocol repeats,
#'   \item average technical replicates,
#'   \item apply a 33-point boxcar filter to selected channels,
#'   \item apply baseline correction to selected channels.
#' }
#'
#' Replicate averaging preserves all channels by default, including the
#' photocell channel, so stimulus alignment and plotting remain available after
#' preprocessing.
#'
#' @param x Either a file path to a `.abf` file, or a `zeus_abf_raw` object.
#' @param ... Passed to `readABF::readABF()` only when `x` is a file path.
#' @param add_stim Logical; if `TRUE`, attach stimulus protocol metadata.
#' @param calib Optional calibration table passed to `nd_to_irradiance_log10()`.
#' @param protocol Character string specifying the protocol. Supported:
#'   `"default"`, `"C1"`, `"C0"`.
#' @param nd_start,nd_end,nd_step,repeats_per_level,n_protocol_repeats Protocol
#'   settings used by `"default"` and `"C1"`. `"C0"` uses a fixed custom
#'   sequence.
#' @param nd_descending Logical; if `TRUE`, protocol ND values are assigned in
#'   descending order.
#' @param treatment_group Treatment label.
#' @param treatment_group_custom Custom label used when
#'   `treatment_group = "user_input"`.
#' @param date_of_fertilization Date of fertilization for the sample.
#' @param erg_age Age at ERG, such as `"Larval"` or `"Adult"`.
#' @param drop_unmatched_sweeps Logical; if `TRUE`, remove sweeps that do not
#'   map to the protocol table.
#' @param average_replicates Logical; if `TRUE`, collapse repeated sweeps into
#'   technical replicate mean waveforms.
#' @param replicate_channel_filter Optional character string specifying which
#'   channel(s) to average. Default is `NULL`, which preserves all channels.
#' @param sweeps_per_replicate Number of repeated sweeps per technical replicate.
#'   Default is `4`.
#' @param add_protocol_repeat Logical; if `TRUE`, add `protocol_repeat` and
#'   `protocol_repeat_rev` columns when protocol metadata are available.
#' @param apply_boxcar Logical; if `TRUE`, apply the validated 33-point boxcar
#'   to channels matching `boxcar_channel_pattern`.
#' @param boxcar_channel_pattern Character string identifying channels to boxcar
#'   filter. Default is `"DAM80"`.
#' @param boxcar_k Integer boxcar width. Default is `33`.
#' @param keep_raw_boxcar Logical; if `TRUE`, preserve the pre-filter signal in
#'   `value_raw`.
#' @param processing_mode Character string controlling replicate-level
#'   correction. One of `"none"`, `"stimresp"`, or `"analysis"`.
#' @param stimresp_offset Optional numeric constant offset subtracted from the
#'   smoothed waveform in `"stimresp"` mode. If `NULL`, no constant offset is
#'   applied automatically.
#' @param baseline_window Numeric length-2 vector specifying the baseline
#'   interval in seconds for `"analysis"` mode. Default is `c(0.300, 0.400)`.
#' @param baseline_channel_pattern Character string identifying channels to
#'   baseline correct. Default is `"DAM80"`.
#'
#' @return A long-format tibble. Depending on options, output may be sweep-level
#'   or replicate-level.
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
    nd_descending = TRUE,
    treatment_group = NULL,
    treatment_group_custom = NULL,
    date_of_fertilization = NA,
    erg_age = NULL,
    drop_unmatched_sweeps = TRUE,
    average_replicates = FALSE,
    replicate_channel_filter = NULL,
    sweeps_per_replicate = 4,
    add_protocol_repeat = TRUE,
    apply_boxcar = FALSE,
    boxcar_channel_pattern = "DAM80",
    boxcar_k = 33,
    keep_raw_boxcar = TRUE,
    processing_mode = c("none", "stimresp", "analysis"),
    stimresp_offset = NULL,
    baseline_window = c(0.300, 0.400),
    baseline_channel_pattern = "DAM80"
) {

  protocol <- match.arg(protocol)
  processing_mode <- match.arg(processing_mode)

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
      n_protocol_repeats = n_protocol_repeats,
      nd_descending = nd_descending,
      treatment_group = treatment_group,
      treatment_group_custom = treatment_group_custom,
      date_of_fertilization = date_of_fertilization,
      erg_age = erg_age
    )

    if (isTRUE(drop_unmatched_sweeps) && "stim_nd" %in% names(df_long)) {
      df_long <- df_long |>
        dplyr::filter(!is.na(.data$stim_nd))
    }

    if (isTRUE(add_protocol_repeat) && "stim_order" %in% names(df_long)) {
      sweeps_per_protocol <- length(
        unique(
          make_protocol_table(
            protocol = protocol,
            nd_start = nd_start,
            nd_end = nd_end,
            nd_step = nd_step,
            repeats_per_level = repeats_per_level,
            n_protocol_repeats = 1
          )$stim_order
        )
      )

      if (is.finite(sweeps_per_protocol) && sweeps_per_protocol > 0) {
        sweep_map <- df_long |>
          dplyr::distinct(.data$sweep, .data$stim_order) |>
          dplyr::arrange(.data$sweep) |>
          dplyr::mutate(
            protocol_repeat = ceiling(.data$stim_order / sweeps_per_protocol)
          )

        max_rep <- max(sweep_map$protocol_repeat, na.rm = TRUE)

        sweep_map <- sweep_map |>
          dplyr::mutate(
            protocol_repeat_rev = max_rep + 1 - .data$protocol_repeat
          )

        df_long <- df_long |>
          dplyr::left_join(sweep_map, by = c("sweep", "stim_order"))
      }
    }
  }

  if (isTRUE(average_replicates)) {
    df_long <- zeus_average_technical_replicates(
      df_long = df_long,
      sweeps_per_replicate = sweeps_per_replicate,
      channel_filter = replicate_channel_filter
    )

    if (isTRUE(apply_boxcar)) {
      df_long <- zeus_boxcar_filter(
        df_long = df_long,
        channel_pattern = boxcar_channel_pattern,
        k = boxcar_k,
        keep_raw = keep_raw_boxcar
      )
    }

    if (processing_mode == "analysis") {
      df_long <- zeus_baseline_correct(
        df_long = df_long,
        baseline_window = baseline_window,
        keep_raw = keep_raw_boxcar,
        channel_pattern = baseline_channel_pattern
      )
    }

    if (processing_mode == "stimresp" && !is.null(stimresp_offset)) {
      if (!is.numeric(stimresp_offset) || length(stimresp_offset) != 1 || is.na(stimresp_offset)) {
        stop("`stimresp_offset` must be a single non-missing numeric value.", call. = FALSE)
      }

      if ("channel" %in% names(df_long)) {
        idx <- grepl(boxcar_channel_pattern, df_long$channel)
      } else {
        idx <- rep(TRUE, nrow(df_long))
      }

      df_long$value[idx] <- df_long$value[idx] - stimresp_offset
    }
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
      sweep = as.integer(s),
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
      sweep = as.integer(s),
      time = time,
      .before = 1
    )
  })
}

# Averaging Technical Replicates ------------------------------------------

#' Average technical replicates into replicate-level waveforms
#'
#' @description
#' Collapses repeated sweeps into technical replicate mean waveforms. By default,
#' all channels are preserved, including photocell channels. If a
#' `channel_filter` is supplied, only matching channels are averaged.
#'
#' @param df_long ZEUS long-format data.
#' @param sweeps_per_replicate Number of sweeps per technical replicate.
#' @param channel_filter Optional channel to average. Default is `NULL`.
#'
#' @return A long-format tibble with one waveform per technical replicate.
#' @keywords internal
zeus_average_technical_replicates <- function(df_long,
                                              sweeps_per_replicate = 4,
                                              channel_filter = NULL) {

  required_cols <- c("sweep", "time", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_work <- df_long

  if (!is.null(channel_filter)) {
    if (!"channel" %in% names(df_work)) {
      stop("`channel_filter` supplied but no `channel` column found.", call. = FALSE)
    }

    df_work <- df_work |>
      dplyr::filter(.data$channel == channel_filter)
  }

  if (nrow(df_work) == 0) {
    stop("No rows remained after channel filtering.", call. = FALSE)
  }

  grouping_keys <- c("stim_nd")
  grouping_keys <- c(grouping_keys, intersect(
    c(
      "stim_irradiance_log10",
      "wavelength",
      "treatment_group",
      "date_of_fertilization",
      "erg_age",
      "protocol_repeat",
      "protocol_repeat_rev"
    ),
    names(df_work)
  ))

  rep_map <- df_work |>
    dplyr::distinct(.data$sweep, dplyr::across(dplyr::all_of(grouping_keys))) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(grouping_keys)), .data$sweep) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) |>
    dplyr::mutate(
      tech_rep = ceiling(dplyr::row_number() / sweeps_per_replicate)
    ) |>
    dplyr::ungroup()

  avg_group_keys <- c(grouping_keys, "tech_rep", "time")
  if ("channel" %in% names(df_work)) avg_group_keys <- c(avg_group_keys, "channel")
  if ("units" %in% names(df_work)) avg_group_keys <- c(avg_group_keys, "units")

  df_work |>
    dplyr::left_join(rep_map, by = c("sweep", grouping_keys)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(avg_group_keys))) |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      sweep = .data$tech_rep
    ) |>
    dplyr::relocate(.data$sweep, .before = .data$time)
}

# Baseline Correction -----------------------------------------------------
#' Apply baseline correction to long-format ERG data
#'
#' @param df_long ZEUS long-format data.
#' @param baseline_window Numeric length-2 vector in seconds.
#' @param keep_raw Logical; if `TRUE`, preserve original signal in `value_raw`.
#'
#' @return Baseline-corrected data.
#' @export
zeus_baseline_correct <- function(df_long,
                                  baseline_window = c(0.300, 0.400),
                                  keep_raw = TRUE,
                                  channel_pattern = "DAM80") {

  required_cols <- c("sweep", "time", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_out <- df_long

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_out)) {
    df_out$value_raw <- df_out$value
  }

  if ("channel" %in% names(df_out)) {
    idx <- grepl(channel_pattern, df_out$channel)
  } else {
    idx <- rep(TRUE, nrow(df_out))
  }

  if (!any(idx)) {
    return(df_out)
  }

  grouping_keys <- intersect(
    c("sweep", "stim_nd", "tech_rep", "wavelength", "protocol_repeat", "protocol_repeat_rev", "channel"),
    names(df_out)
  )

  df_target <- df_out[idx, , drop = FALSE]
  df_other  <- df_out[!idx, , drop = FALSE]

  df_target <- df_target |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) |>
    dplyr::mutate(
      baseline = mean(
        .data$value[
          .data$time >= baseline_window[1] &
            .data$time <= baseline_window[2]
        ],
        na.rm = TRUE
      ),
      value = .data$value - .data$baseline
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$baseline)

  dplyr::bind_rows(df_target, df_other) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(grouping_keys)), .data$time)
}

# Boxcar filter -----------------------------------------------------------

#' Apply a legacy-style boxcar filter to long-format ERG data
#'
#' @description
#' Applies a fixed-width running-average (boxcar) filter to the `value` column
#' within each replicate waveform for channels matching `channel_pattern`.
#'
#' @param df_long ZEUS long-format data.
#' @param channel_pattern Character string identifying channel(s) to filter.
#' @param k Integer window size. Default is `33`.
#' @param keep_raw Logical; if `TRUE`, preserve original signal in `value_raw`.
#'
#' @return Filtered data.
#' @export
zeus_boxcar_filter <- function(df_long,
                               channel_pattern = "DAM80",
                               k = 33,
                               keep_raw = TRUE) {

  required_cols <- c("sweep", "time", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_out <- df_long

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_out)) {
    df_out$value_raw <- df_out$value
  }

  if ("channel" %in% names(df_out)) {
    idx <- grepl(channel_pattern, df_out$channel)
  } else {
    idx <- rep(TRUE, nrow(df_out))
  }

  if (!any(idx)) {
    return(df_out)
  }

  legacy_boxcar_fixed <- function(x, k = 33) {
    n <- length(x)
    out <- numeric(n)
    half_k <- (k - 1L) %/% 2L

    for (i in seq_len(n)) {
      start_idx <- i - half_k
      end_idx <- i + half_k

      if (start_idx < 1L) {
        start_idx <- 1L
        end_idx <- min(k, n)
      }

      if (end_idx > n) {
        end_idx <- n
        start_idx <- max(1L, n - k + 1L)
      }

      out[i] <- mean(x[start_idx:end_idx], na.rm = TRUE)
    }

    out
  }

  grouping_keys <- intersect(
    c("sweep", "stim_nd", "tech_rep", "wavelength", "protocol_repeat", "protocol_repeat_rev", "channel"),
    names(df_out)
  )

  df_target <- df_out[idx, , drop = FALSE]
  df_other  <- df_out[!idx, , drop = FALSE]

  df_target <- df_target |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::mutate(
      value = legacy_boxcar_fixed(.data$value, k = k)
    ) |>
    dplyr::ungroup()

  dplyr::bind_rows(df_target, df_other) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(grouping_keys)), .data$time)
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
#' @param nd_descending Logical; if `TRUE`, order ND levels from high to low.
#'
#' @return A data.frame describing the sweep protocol.
#' @keywords internal
make_protocol_default <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10,
    nd_descending = TRUE
) {
  nd_levels <- seq(nd_start, nd_end, by = nd_step)

  if (isTRUE(nd_descending)) {
    nd_levels <- rev(nd_levels)
  }

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
#' @param nd_descending Logical; if `TRUE`, order ND levels from high to low.
#'
#' @return A data.frame describing the sweep protocol.
#' @keywords internal
make_protocol_C1 <- function(
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10,
    nd_descending = TRUE
) {
  nd_levels <- seq(nd_start, nd_end, by = nd_step)

  if (isTRUE(nd_descending)) {
    nd_levels <- rev(nd_levels)
  }

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
#' @param nd_descending Logical; if `TRUE`, order ND levels from high to low
#'   for default/C1 protocols.
#'
#' @return A protocol table with `stim_order`, `stim_type`, `stim_nd`, and `wavelength`.
#' @keywords internal
make_protocol_table <- function(
    protocol = c("default", "C1", "C0"),
    nd_start = 3.0,
    nd_end = 6.0,
    nd_step = 0.5,
    repeats_per_level = 4,
    n_protocol_repeats = 10,
    nd_descending = TRUE
) {
  protocol <- match.arg(protocol)

  switch(
    protocol,
    default = make_protocol_default(
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats,
      nd_descending = nd_descending
    ),
    C1 = make_protocol_C1(
      nd_start = nd_start,
      nd_end = nd_end,
      nd_step = nd_step,
      repeats_per_level = repeats_per_level,
      n_protocol_repeats = n_protocol_repeats,
      nd_descending = nd_descending
    ),
    C0 = make_protocol_C0(
      repeats_per_level = repeats_per_level
    )
  )
}


# Allowed Treatment values -----------------------------------------------------

#' Allowed treatment group values
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
#' @return A character vector of allowed ERG age values.
#' @keywords internal
.allowed_erg_age <- function() {
  c("Larval", "Adult")
}

# Treatment and Age Metadata Helper --------------------------------------------

#' Construct sample-level metadata for ERG recordings
#'
#' @param treatment_group Character string specifying treatment group.
#' @param treatment_group_custom Character string used when `treatment_group`
#'   = `"user_input"` to specify a custom treatment label.
#' @param date_of_fertilization Date of fertilization for the fish.
#' @param erg_age Character string describing developmental stage at ERG.
#'
#' @return A data frame containing validated sample metadata.
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
#' @param nd_descending Logical; if `TRUE`, order ND levels from high to low
#'   for default/C1 protocols.
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
    nd_descending = TRUE,
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
    n_protocol_repeats = n_protocol_repeats,
    nd_descending = nd_descending
  )

  sweeps_present <- sort(unique(as.integer(df_long$sweep)))
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
    dplyr::mutate(sweep = as.integer(.data$sweep)) |>
    dplyr::left_join(
      stim_tbl,
      by = "sweep",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      treatment_group = sample_meta$treatment_group,
      date_of_fertilization = sample_meta$date_of_fertilization,
      erg_age = sample_meta$erg_age
    ) |>
    dplyr::filter(!is.na(.data$stim_nd))
}

