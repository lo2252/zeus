#' Package-level imports and NSE global variable declarations
#'
#' Internal package annotations used to keep the generated namespace in sync and
#' to silence `R CMD check` notes for non-standard evaluation in dplyr-powered
#' code paths.
#'
#' @keywords internal
#' @noRd
#' @importFrom rlang .data
#' @importFrom stats ave
NULL

utils::globalVariables(c(
  "baseline",
  "grouping_value",
  "grouping_variable",
  "keep",
  "max_peak_mv",
  "mean_latency_ms",
  "mean_peak_mv",
  "median_peak_mv",
  "min_peak_mv",
  "n",
  "peak_type",
  "protocol_id",
  "sd_latency_ms",
  "sd_peak_mv",
  "sem_peak_mv",
  "signal",
  "stim_index",
  "stim_nd",
  "summary_type",
  "value",
  "wavelength",
  "wavelength_label"
))
