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

#' Import and preprocess ABF electrophysiology data into ZEUS format
#'
#' @description
#' Imports Axon Binary Format (ABF) electrophysiology data and converts it into
#' a standardized long-format data frame suitable for downstream analysis within
#' the ZEUS framework.
#'
#' The function supports optional preprocessing steps including baseline
#' correction and application of a legacy-style boxcar filter, as well as
#' automatic assignment of stimulus protocol information and experimental
#' metadata.
#'
#' @details
#' The preprocessing pipeline proceeds in the following order:
#' \enumerate{
#'   \item Import ABF file using \code{read_abf_raw()}.
#'   \item Convert to long-format data using \code{abf_as_df_long()}.
#'   \item (Optional) Apply baseline correction using \code{zeus_baseline_correct()}.
#'   \item (Optional) Apply a legacy 33-point boxcar filter using \code{zeus_boxcar_filter()}.
#'   \item (Optional) Append stimulus protocol and metadata using
#'         \code{add_stimulus_cols_protocol()}.
#' }
#'
#' The boxcar filter is implemented as a fixed-width 33-point running average,
#' corresponding to a 16.5 ms window when the sampling interval is 0.5 ms.
#' Edge behavior is preserved to match legacy electrophysiology workflows,
#' resulting in a constant-valued region at the beginning and end of each trace.
#'
#' Raw (unfiltered) signal values can optionally be preserved in a separate
#' column (\code{value_raw}) when filtering is applied.
#'
#' @param x Either a file path to an ABF file or a \code{zeus_abf_raw} object.
#'
#' @param ... Additional arguments passed to \code{read_abf_raw()}.
#'
#' @param add_stim Logical; if \code{TRUE}, appends stimulus protocol
#'   information to the output. Default is \code{TRUE}.
#'
#' @param calib Optional calibration data used for mapping stimulus ND values
#'   to irradiance. If \code{NULL}, a default calibration is used.
#'
#' @param protocol Character string specifying the stimulus protocol.
#'   Options include:
#'   \itemize{
#'     \item \code{"default"}: ND sweep from \code{nd_start} to \code{nd_end}
#'     \item \code{"C1"}: constant wavelength protocol
#'     \item \code{"C0"}: spectral protocol with wavelength-specific ND ordering
#'   }
#'
#' @param nd_start Starting neutral density (ND) value. Default is \code{3.0}.
#'
#' @param nd_end Ending ND value. Default is \code{6.0}.
#'
#' @param nd_step Step size between ND levels. Default is \code{0.5}.
#'
#' @param repeats_per_level Number of sweeps per ND level. Default is \code{4}.
#'
#' @param n_protocol_repeats Number of full protocol repetitions. Default is \code{10}.
#'
#' @param nd_descending Logical; if \code{TRUE}, ND levels are ordered from high
#'   to low intensity. Default is \code{TRUE}.
#'
#' @param treatment_group Character string specifying treatment condition.
#'
#' @param treatment_group_custom Optional custom label when
#'   \code{treatment_group = "user_input"}.
#'
#' @param date_of_fertilization Date of fertilization for the sample. Can be
#'   coerced using \code{as.Date()} or parsed externally.
#'
#' @param erg_age Character string specifying developmental stage (e.g.,
#'   \code{"Larval"}, \code{"Adult"}).
#'
#' @param apply_baseline Logical; if \code{TRUE}, performs baseline correction
#'   prior to filtering. Default is \code{FALSE}.
#'
#' @param baseline_window Numeric vector of length 2 specifying the time window
#'   (in seconds) used for baseline estimation. Default is \code{c(0, 0.3)}.
#'
#' @param apply_boxcar Logical; if \code{TRUE}, applies a boxcar filter to the
#'   signal. Default is \code{TRUE}.
#'
#' @param boxcar_channel_pattern Character string used to identify channels to
#'   filter (e.g., \code{"DAM80"}). Default is \code{"DAM80"}.
#'
#' @param boxcar_k Integer window size for the running average. Default is
#'   \code{33}. If \code{NULL}, defaults to 33.
#'
#' @param keep_raw_boxcar Logical; if \code{TRUE}, preserves original signal in
#'   \code{value_raw}. Default is \code{TRUE}.
#'
#' @return A tibble in ZEUS long format containing:
#' \itemize{
#'   \item \code{sweep}: sweep index
#'   \item \code{time}: time in seconds
#'   \item \code{channel}: channel name
#'   \item \code{value}: processed signal (baseline-corrected and/or filtered)
#'   \item \code{value_raw}: original signal (if \code{keep_raw_boxcar = TRUE})
#'   \item stimulus protocol columns (e.g., \code{stim_nd},
#'         \code{stim_irradiance_log10})
#'   \item metadata columns (e.g., \code{treatment_group}, \code{erg_age})
#' }
#'
#' @examples
#' \dontrun{
#' df <- zeus_import(
#'   "example.abf",
#'   protocol = "C1",
#'   treatment_group = "SYS Water",
#'   erg_age = "Larval"
#' )
#' }
#'
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
    apply_baseline = TRUE,
    baseline_window = c(0, 0.3),
    apply_boxcar = TRUE,
    boxcar_channel_pattern = "DAM80",
    boxcar_k = NULL,
    keep_raw_boxcar = TRUE
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

  # Apply baseline
  if (isTRUE(apply_baseline)) {
    df_long <- zeus_baseline_correct(
      df_long = df_long,
      baseline_window = baseline_window
    )
  }

  # Apply Boxcar
  if (isTRUE(apply_boxcar)) {
    df_long <- zeus_boxcar_filter(
      df_long,
      channel_pattern = boxcar_channel_pattern,
      k = if (is.null(boxcar_k)) 33 else boxcar_k,
      keep_raw = keep_raw_boxcar
    )
  }

  # Add stim protocol
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


# Baseline Correction -----------------------------------------------------
#' Apply baseline correction to long-format ERG data
#'
#' @description
#' Performs sweep-wise baseline correction by subtracting the mean signal value
#' within a user-defined baseline window from all time points in each sweep and
#' channel.
#'
#' This function is intended for use on ZEUS long-format data and is typically
#' applied before smoothing and feature extraction.
#'
#' @details
#' Baseline correction is performed independently for each combination of
#' `sweep` and `channel`.
#'
#' For each sweep-channel pair, the baseline is defined as the mean of `value`
#' within the interval specified by `baseline_window`. That baseline value is
#' then subtracted from the full trace.
#'
#' If `keep_raw = TRUE`, the original uncorrected signal is preserved in a
#' separate column named `value_raw`, while the corrected signal is stored in
#' `value`.
#'
#' This function does not smooth the signal and does not modify stimulus or
#' metadata columns.
#'
#' @param df_long A ZEUS long-format data frame. Must contain at minimum
#'   `sweep`, `time`, `channel`, and `value`.
#'
#' @param baseline_window Numeric vector of length 2 giving the start and end
#'   of the baseline window in the same units as `df_long$time`. Default is
#'   `c(0, 0.3)`.
#'
#' @param keep_raw Logical; if `TRUE`, preserves the original uncorrected
#'   signal in a separate column named `value_raw`. Default is `TRUE`.
#'
#' @return A data frame with baseline-corrected values stored in `value`.
#'   If `keep_raw = TRUE`, the original signal is also retained in `value_raw`.
#' @export
zeus_baseline_correct <- function(df_long,
                                  baseline_window = c(0, 0.3),
                                  keep_raw = TRUE) {

  required_cols <- c("sweep", "time", "channel", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(baseline_window) ||
      length(baseline_window) != 2 ||
      any(is.na(baseline_window)) ||
      baseline_window[1] >= baseline_window[2]) {
    stop(
      "`baseline_window` must be a numeric vector of length 2 with start < end.",
      call. = FALSE
    )
  }

  df_out <- df_long

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_out)) {
    df_out$value_raw <- df_out$value
  }

  df_out |>
    dplyr::group_by(.data$sweep, .data$channel) |>
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
}

# Boxcar filter -----------------------------------------------------------

#' Apply a legacy-style boxcar filter to long-format ERG data
#'
#' @description
#' Applies a fixed-width running-average (boxcar) filter to the `value` column
#' within each sweep for channels matching `channel_pattern`.
#'
#' This function is designed to reproduce a legacy electrophysiology filtering
#' workflow using a 33-point running average. When the sampling interval is
#' 0.5 ms, this corresponds to a 16.5 ms smoothing window.
#'
#' @details
#' The filter is applied independently within each sweep and channel.
#'
#' To preserve compatibility with the legacy workflow, edge handling uses a
#' fixed-width window: near the beginning of a trace, the first full window is
#' reused until the averaging window can slide normally; near the end of a
#' trace, the last full window is reused similarly. This may produce a flat
#' region at the start and end of the filtered waveform and is intentional.
#'
#' If `keep_raw = TRUE`, the original unsmoothed signal is preserved in a
#' separate column named `value_raw`, while the filtered signal is stored in
#' `value`.
#'
#' @param df_long A ZEUS long-format data frame. Must contain at minimum
#'   `sweep`, `time`, `channel`, and `value`.
#'
#' @param channel_pattern Character string used to identify the channel(s) to
#'   filter. Default is `"DAM80"`.
#'
#' @param k Integer window size for the running average. Default is `33`.
#'   Must be an odd positive integer.
#'
#' @param keep_raw Logical; if `TRUE`, the original signal is preserved in
#'   `value_raw`. Default is `TRUE`.
#'
#' @return A data frame with the same structure as `df_long`, with filtered
#'   values stored in `value`. If `keep_raw = TRUE`, the original signal is
#'   retained in `value_raw`.
#' @keywords internal
zeus_boxcar_filter <- function(df_long,
                               channel_pattern = "DAM80",
                               k = 33,
                               keep_raw = TRUE) {

  required_cols <- c("sweep", "time", "channel", "value")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1) {
    stop("`k` must be a single integer >= 1.", call. = FALSE)
  }

  k <- as.integer(k)

  if (k %% 2 == 0) {
    stop("`k` must be odd for the legacy symmetric running average.", call. = FALSE)
  }

  idx <- grepl(channel_pattern, df_long$channel)

  if (!any(idx)) {
    if (isTRUE(keep_raw) && !"value_raw" %in% names(df_long)) {
      df_long$value_raw <- df_long$value
    }
    return(df_long)
  }

  df_target <- df_long[idx, , drop = FALSE]
  df_other  <- df_long[!idx, , drop = FALSE]

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_target)) {
    df_target$value_raw <- df_target$value
  }

  if (isTRUE(keep_raw) && !"value_raw" %in% names(df_other)) {
    df_other$value_raw <- df_other$value
  }

  legacy_boxcar <- function(x, k) {
    n <- length(x)
    out <- numeric(n)
    half_k <- (k - 1L) %/% 2L

    for (i in seq_len(n)) {
      start_idx <- i - half_k
      end_idx   <- i + half_k

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

  split_groups <- split(
    df_target,
    interaction(df_target$sweep, df_target$channel, drop = TRUE)
  )

  filtered_groups <- lapply(split_groups, function(dat) {
    dat <- dat[order(dat$time), , drop = FALSE]
    dat$value <- legacy_boxcar(dat$value, k = k)
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

