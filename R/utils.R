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
#' Smoothing helper used for ZEUS traces.
#' If `n = 1`, the input is returned unchanged.
#' If `n > 1`, a centered moving average is applied with endpoint extension.
#'
#' @param x Numeric vector.
#' @param n Integer smoothing window size.
#'
#' @return Numeric vector of the same length as the input.
#' @keywords internal
run_mean <- function(x, n = 1L) {
  n <- as.integer(n)

  if (n <= 1L) {
    return(x)
  }

  zoo::rollmean(x, k = n, fill = "extend")
}