# Time Conversion Function ------------------------------------------------

#' Convert time values to milliseconds when needed
#'
#' Internal helper to standardize time units across ZEUS functions.
#' This helper uses a simple rule:
#' - If the maximum time is small (<= 20), treat the input as seconds
#'   and convert to milliseconds.
#' - Otherwise, assume the input is already in milliseconds.
#'
#' @param x Numeric vector of time values.
#'
#' @return Numeric vector in milliseconds.
#' @keywords internal
zeus_time_to_ms <- function(x) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }

  xmax <- suppressWarnings(max(x, na.rm = TRUE))

  if (is.finite(xmax) && xmax <= 20) {
    x * 1000
  } else {
    x
  }
}


# Stimulus Labels ---------------------------------------------------------

#' Build StimResp labels
#'
#' Creates the labels used for averaged StimResp traces.
#'
#' @param wavelength Character or numeric wavelength value.
#' @param stim_nd Numeric ND value.
#' @param stim_index Optional integer stimulus index.
#' @param protocol_id Optional protocol id, such as `"C0"` or `"C1"`.
#'
#' @return Character scalar label.
#' @keywords internal
make_stimresp_label <- function(wavelength,
                                stim_nd,
                                stim_index = NULL,
                                protocol_id = NULL) {
  # Convert wavelength to character so both numeric wavelengths
  # and labels like "White" are handled consistently.
  wl <- as.character(wavelength)

  # Format ND values with one decimal place.
  nd <- formatC(stim_nd, format = "f", digits = 1)

  # Special C0 labeling for the two 650 nm blocks.
  if (identical(protocol_id, "C0")) {
    suffix <- ""

    # First 650 block gets "A"
    if (!is.null(stim_index) && stim_index >= 1L && stim_index <= 7L) {
      suffix <- "A"
    }

    # Second 650 block gets "B"
    if (!is.null(stim_index) && stim_index >= 36L && stim_index <= 42L) {
      suffix <- "B"
    }

    paste0(wl, suffix, " ", nd)
  } else {
    # Other protocols use the simpler label format.
    paste0(wl, " ", nd)
  }
}


# Centered Running Mean ---------------------------------------------------

#' Compute a centered running mean
#'
#' Smoothing helper used for StimResp traces.
#' If `n = 1`, the input is returned unchanged.
#' If `n > 1`, a centered moving average is applied.
#'
#' @param x Numeric vector.
#' @param n Integer smoothing window size.
#'
#' @return Numeric vector of approximately the same length as the input.
#' @keywords internal
run_mean <- function(x, n = 1L) {
  # Convert window size to integer.
  n <- as.integer(n)

  # No smoothing.
  if (n <= 1L) {
    return(x)
  }

  # Apply a centered moving average.
  stats::filter(x, rep(1 / n, n), sides = 2) |>
    as.numeric()
}


# Validation of long-form ABF data ----------------------------------------

#' Validate the long-form structure of ABF data
#'
#' Confirms that `abf_as_df_long()` returned the minimum columns ZEUS expects.
#'
#' Required columns:
#' - `sweep`
#' - `time`
#' - `channel`
#' - `value`
#'
#' This function does not modify the data. It only checks structure.
#'
#' @param df_long Data frame returned by `abf_as_df_long()`.
#'
#' @return Invisibly returns `df_long`.
#' @keywords internal
validate_long_abf_df <- function(df_long) {
  needed <- c("sweep", "time", "channel", "value")
  missing_cols <- setdiff(needed, names(df_long))

  if (length(missing_cols) > 0L) {
    stop(
      "`abf_as_df_long()` must return columns: ",
      paste(needed, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(df_long)
}


# Individual Channel Extraction -------------------------------------------

#' Extract one channel from long-form ABF data
#'
#' Filters the long ABF data to a single channel, such as:
#' - `"ERG DAM80"`
#' - `"Photocell"`
#' - numeric channel ids
#'
#' Used to separate ERG traces from photocell traces before building
#' StimResp outputs.
#'
#' @param df_long Long-form ABF data.
#' @param channel_value Channel identifier to keep.
#'
#' @return Filtered data frame containing only rows from the requested channel.
#' @keywords internal
extract_channel_trace_df <- function(df_long, channel_value) {
  out <- df_long |>
    dplyr::filter(.data$channel == channel_value)

  if (nrow(out) == 0L) {
    stop("No rows found for channel '", channel_value, "'.", call. = FALSE)
  }

  out
}
