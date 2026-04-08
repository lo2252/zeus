# Importing Raw ERG -------------------------------------------------------

#' Import a raw ABF file into a ZEUS raw object
#'
#' @description
#' Reads an Axon Binary Format (`.abf`) file using `readABF::readABF()` and
#' wraps the result in a `zeus_abf_raw` object for downstream ZEUS processing.
#'
#' @param path Path to a `.abf` file.
#' @param ... Additional arguments passed to `readABF::readABF()`.
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

#' Read ABF data into ZEUS
#'
#' @description
#' User-facing import wrapper for ZEUS ABF files.
#'
#' This function can:
#' \enumerate{
#'   \item return the raw imported ABF object, or
#'   \item build a StimResp-style processed object using a supported protocol.
#' }
#'
#' Currently supported protocols:
#' \itemize{
#'   \item `"raw"`: return only the raw imported ABF object
#'   \item `"C0"`: spectral protocol
#'   \item `"C1"`: white-light protocol
#' }
#'
#' @param x File path to a `.abf` file or a `zeus_abf_raw` object.
#' @param protocol One of `"raw"`, `"C0"`, or `"C1"`.
#' @param erg_channel ERG channel identifier.
#' @param pc_channel Photocell channel identifier.
#' @param repeats_per_stim Integer number of technical replicates per stimulus.
#' @param expected_stim Integer number of stimulus conditions.
#' @param exclude_noisy Logical.
#' @param noise_threshold Numeric threshold for noisy-trace exclusion.
#' @param zero_baseline Logical.
#' @param baseline_window_ms Numeric length-2 vector in milliseconds.
#' @param smooth_n Integer running-mean window size.
#' @param align_to_stimulus One of `"protocol"` or `"photocell"`.
#' @param photocell_baseline_window_ms Numeric length-2 vector used for photocell baseline estimation.
#' @param photocell_threshold_frac Fraction of photocell rise used for onset detection.
#' @param ... Additional arguments passed only to `read_abf_raw()` when `x` is a file path.
#'
#' @return A `zeus_abf_raw` object if `protocol = "raw"`, otherwise a
#'   `zeus_stimresp` object.
#' @export
zeus_read_abf <- function(x,
                          protocol = c("raw", "C0", "C1"),
                          erg_channel = "ERG DAM80",
                          pc_channel = "Photocell",
                          repeats_per_stim = 4L,
                          expected_stim = 70L,
                          exclude_noisy = FALSE,
                          noise_threshold = 1.5,
                          zero_baseline = TRUE,
                          baseline_window_ms = c(300, 400),
                          smooth_n = 1L,
                          align_to_stimulus = c("protocol", "photocell"),
                          photocell_baseline_window_ms = c(0, 300),
                          photocell_threshold_frac = 0.5,
                          ...) {
  protocol <- match.arg(protocol)
  align_to_stimulus <- match.arg(align_to_stimulus)

  repeats_per_stim <- as.integer(repeats_per_stim)
  expected_stim <- as.integer(expected_stim)
  smooth_n <- as.integer(smooth_n)

  if (repeats_per_stim < 1L) {
    stop("'repeats_per_stim' must be >= 1.", call. = FALSE)
  }

  if (expected_stim < 1L) {
    stop("'expected_stim' must be >= 1.", call. = FALSE)
  }

  if (!is.logical(exclude_noisy) || length(exclude_noisy) != 1L || is.na(exclude_noisy)) {
    stop("'exclude_noisy' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(noise_threshold) || length(noise_threshold) != 1L || is.na(noise_threshold)) {
    stop("'noise_threshold' must be a single numeric value.", call. = FALSE)
  }

  if (!is.logical(zero_baseline) || length(zero_baseline) != 1L || is.na(zero_baseline)) {
    stop("'zero_baseline' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(baseline_window_ms) || length(baseline_window_ms) != 2L) {
    stop("'baseline_window_ms' must be a numeric vector of length 2.", call. = FALSE)
  }

  if (smooth_n < 1L) {
    stop("'smooth_n' must be >= 1.", call. = FALSE)
  }

  if (!is.numeric(photocell_baseline_window_ms) || length(photocell_baseline_window_ms) != 2L) {
    stop(
      "'photocell_baseline_window_ms' must be a numeric vector of length 2.",
      call. = FALSE
    )
  }

  if (!is.numeric(photocell_threshold_frac) ||
      length(photocell_threshold_frac) != 1L ||
      is.na(photocell_threshold_frac) ||
      photocell_threshold_frac <= 0) {
    stop(
      "'photocell_threshold_frac' must be a single positive numeric value.",
      call. = FALSE
    )
  }

  raw_obj <- if (inherits(x, "zeus_abf_raw")) {
    x
  } else if (is.character(x) && length(x) == 1L) {
    read_abf_raw(path = x, ...)
  } else {
    stop("'x' must be a file path or a 'zeus_abf_raw' object.", call. = FALSE)
  }

  if (identical(protocol, "raw")) {
    return(raw_obj)
  }

  build_stimresp(
    x = raw_obj,
    protocol = protocol,
    erg_channel = erg_channel,
    pc_channel = pc_channel,
    repeats_per_stim = repeats_per_stim,
    expected_stim = expected_stim,
    exclude_noisy = exclude_noisy,
    noise_threshold = noise_threshold,
    zero_baseline = zero_baseline,
    baseline_window_ms = baseline_window_ms,
    smooth_n = smooth_n,
    align_to_stimulus = align_to_stimulus,
    photocell_baseline_window_ms = photocell_baseline_window_ms,
    photocell_threshold_frac = photocell_threshold_frac
  )
}

# Convert to long DF ------------------------------------------------------

#' Convert readABF output to a standardized long data frame
#'
#' @param raw_abf An object of class `.abf`.
#'
#' @return A tibble in long format with columns:
#'   `sweep`, `time`, `channel`, `units`, and `value`.
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

#' Convert readABF output to a standardized wide data frame
#'
#' @param raw_abf An object of class `ABF`.
#'
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
