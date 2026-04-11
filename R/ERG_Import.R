# Importing Raw ERG -------------------------------------------------------

#' Import a raw ABF file into a ZEUS raw object
#'
#' @description
#' Reads an Axon Binary Format (`.abf`) file using [readABF::readABF()] and
#' wraps the result in a `zeus_abf_raw` object for downstream ZEUS processing.
#'
#' @param path Path to a `.abf` file.
#' @param ... Additional arguments passed to [readABF::readABF()].
#'
#' @return An object of class `zeus_abf_raw`.
#' @export
read_abf_raw <- function(path, ...) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("'path' must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

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

# Read ABF wrapper --------------------------------------------------------

#' Read an ABF file and build a ZEUS stimulus-response object
#'
#' @description
#' Wrapper around [build_stimresp()] that reads an ABF file, extracts ERG and
#' optional photocell channels, labels stimulus structure, and returns a
#' structured `zeus_stimresp` object containing raw, replicate-level, and
#' averaged stimulus-response outputs.
#'
#' @param path File path to an ABF file.
#' @param protocol Protocol id. One of `"C0"` or `"C1"`.
#' @param erg_channel Name of the ERG channel in the ABF file. Default is
#'   `"ERG DAM80"`.
#' @param pc_channel Name of the photocell channel in the ABF file. Default is
#'   `"Photocell"`.
#' @param zero_baseline Logical; apply baseline subtraction at the
#'   technical-replicate level before averaging. Default is `FALSE`.
#' @param smooth_n Integer running-mean window size applied to replicate-level
#'   traces before any downstream processing. Default is `1L`.
#' @param align_to_stimulus One of `"protocol"` or `"photocell"`.
#' @param stimresp_zero_baseline Logical; apply baseline subtraction to
#'   replicate-level traces before averaging into stimulus-response traces.
#'   Default is `TRUE`.
#' @param stimresp_baseline_window_ms Numeric length-2 vector for baseline
#'   subtraction. Default is `c(300, 400)`.
#' @param stimresp_exclude_noisy Logical; exclude noisy technical replicates
#'   before averaging. Default is `TRUE`.
#' @param stimresp_noise_window_ms Numeric length-2 vector for replicate-level
#'   noise estimation. Default is `c(300, 1000)`.
#' @param stimresp_noise_ratio_cutoff Numeric SD-ratio threshold used to remove
#'   noisy technical replicates. Default is `1.5`.
#' @param stimresp_runmean_k Integer running-average window size for averaged
#'   stimulus-response traces. Use `1L` for none. Default is `16L`.
#' @param stimresp_drift_correct Logical; stored in the output object for
#'   downstream workflows. Default is `FALSE`.
#' @param stimresp_drift_tail_n Integer number of final points used for
#'   downstream drift estimation. Default is `100L`.
#' @param stimresp_sg_smooth Logical; apply Savitzky-Golay smoothing to averaged
#'   stimulus-response traces. Default is `TRUE`.
#' @param stimresp_sg_n Integer Savitzky-Golay window size for averaged traces.
#'   Default is `101L`.
#' @param stimresp_sg_p Integer Savitzky-Golay polynomial order. Default is `2L`.
#' @param ... Additional arguments passed to [read_abf_raw()].
#'
#' @return An object of class `zeus_stimresp`.
#' @export
zeus_read_abf <- function(path,
                          protocol = c("C0", "C1"),
                          erg_channel = "ERG DAM80",
                          pc_channel = "Photocell",
                          zero_baseline = FALSE,
                          smooth_n = 1L,
                          align_to_stimulus = c("protocol", "photocell"),
                          stimresp_zero_baseline = TRUE,
                          stimresp_baseline_window_ms = c(300, 400),
                          stimresp_exclude_noisy = TRUE,
                          stimresp_noise_window_ms = c(300, 1000),
                          stimresp_noise_ratio_cutoff = 1.5,
                          stimresp_runmean_k = 16L,
                          stimresp_drift_correct = FALSE,
                          stimresp_drift_tail_n = 100L,
                          stimresp_sg_smooth = TRUE,
                          stimresp_sg_n = 101L,
                          stimresp_sg_p = 2L,
                          ...) {
  protocol <- match.arg(protocol)
  align_to_stimulus <- match.arg(align_to_stimulus)

  raw_obj <- read_abf_raw(path, ...)
  raw_abf <- raw_obj$raw

  abf_long <- abf_as_df_long(raw_abf)

  erg_df <- abf_long |>
    dplyr::filter(.data$channel == erg_channel)

  if (nrow(erg_df) == 0L) {
    stop("No ERG data found for channel: ", erg_channel, call. = FALSE)
  }

  pc_df <- abf_long |>
    dplyr::filter(.data$channel == pc_channel)

  if (isTRUE(zero_baseline)) {
    erg_df <- zero_baseline_traces(
      erg_df = erg_df,
      baseline_window_ms = c(300, 400)
    )
  }

  smooth_n <- as.integer(smooth_n)
  if (smooth_n > 1L) {
    erg_df <- smooth_stimresp_traces(
      erg_df = erg_df,
      n = smooth_n
    )
  }

  traces_280 <- label_stimulus_protocol(
    erg_df = erg_df,
    protocol = protocol
  ) |>
    dplyr::mutate(
      time_ms = zeus_time_to_ms(.data$time)
    )

  if (identical(align_to_stimulus, "photocell")) {
    if (nrow(pc_df) == 0L) {
      warning(
        "Photocell alignment requested but no photocell channel was found. Falling back to protocol alignment.",
        call. = FALSE
      )
    } else {
      onset_df <- compute_photocell_onsets(pc_df)

      traces_280 <- add_relative_time_cols(
        df = traces_280,
        onset_df = onset_df,
        by = "sweep",
        onset_col = "stim_onset_ms"
      )
    }
  }

  stimresp_obj <- build_stimresp(
    traces_280 = traces_280,
    stimresp_zero_baseline = stimresp_zero_baseline,
    stimresp_baseline_window_ms = stimresp_baseline_window_ms,
    stimresp_exclude_noisy = stimresp_exclude_noisy,
    stimresp_noise_window_ms = stimresp_noise_window_ms,
    stimresp_noise_ratio_cutoff = stimresp_noise_ratio_cutoff,
    stimresp_runmean_k = stimresp_runmean_k,
    stimresp_sg_smooth = stimresp_sg_smooth,
    stimresp_sg_n = stimresp_sg_n,
    stimresp_sg_p = stimresp_sg_p
  )

  traces_70 <- stimresp_obj$traces_70
  stimresp_qc <- stimresp_obj$stimresp_qc

  out <- list(
    raw = raw_abf,
    path = raw_obj$path,
    protocol = protocol,
    erg_channel = erg_channel,
    pc_channel = pc_channel,
    align_to_stimulus = align_to_stimulus,
    traces_280 = traces_280,
    traces_70 = traces_70,
    stimresp_qc = stimresp_qc,
    stimresp_settings = list(
      stimresp_zero_baseline = stimresp_zero_baseline,
      stimresp_baseline_window_ms = stimresp_baseline_window_ms,
      stimresp_exclude_noisy = stimresp_exclude_noisy,
      stimresp_noise_window_ms = stimresp_noise_window_ms,
      stimresp_noise_ratio_cutoff = stimresp_noise_ratio_cutoff,
      stimresp_runmean_k = stimresp_runmean_k,
      stimresp_drift_correct = stimresp_drift_correct,
      stimresp_drift_tail_n = stimresp_drift_tail_n,
      stimresp_sg_smooth = stimresp_sg_smooth,
      stimresp_sg_n = stimresp_sg_n,
      stimresp_sg_p = stimresp_sg_p
    )
  )

  class(out) <- c("zeus_stimresp", "zeus_abf", class(out))
  out
}

# Convert to long DF ------------------------------------------------------

#' Convert ABF output to a standardized long data frame
#'
#' @param raw_abf An `ABF` object or a `zeus_abf_raw` object.
#'
#' @return A tibble in long format with columns `sweep`, `time`, `channel`,
#'   `units`, and `value`.
#' @keywords internal
abf_as_df_long <- function(raw_abf) {
  if (inherits(raw_abf, "zeus_abf_raw")) {
    raw_abf <- raw_abf$raw
  }

  if (!inherits(raw_abf, "ABF")) {
    stop("`raw_abf` must be an ABF object or a `zeus_abf_raw` object.", call. = FALSE)
  }

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

# Convert to wide DF ------------------------------------------------------

#' Convert ABF output to a standardized wide data frame
#'
#' @param raw_abf An `ABF` object or a `zeus_abf_raw` object.
#'
#' @return A tibble in wide format with one row per time point per sweep.
#' @keywords internal
abf_as_df_wide <- function(raw_abf) {
  if (inherits(raw_abf, "zeus_abf_raw")) {
    raw_abf <- raw_abf$raw
  }

  if (!inherits(raw_abf, "ABF")) {
    stop("`raw_abf` must be an ABF object or a `zeus_abf_raw` object.", call. = FALSE)
  }

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