# Setting Environment -----------------------------------------------------

.zeus_irrad_registry <- new.env(parent = emptyenv())


# Irradiance Calibration Table --------------------------------------------

#' Register irradiance calibration table
#'
#' @param name Calibration id.
#' @param calib Data frame with columns `wavelength`, `stim_nd`, `log_hv`.
#'
#' @return Invisibly returns `calib`.
#' @export
register_zeus_irrad_calibration <- function(name, calib) {
  required <- c("wavelength", "stim_nd", "log_hv")
  missing_cols <- setdiff(required, names(calib))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Calibration table must contain: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  
  calib <- tibble::as_tibble(calib) |>
    dplyr::mutate(wavelength = as.character(.data$wavelength))
  
  assign(name, calib, envir = .zeus_irrad_registry)
  invisible(calib)
}


# Registered Irradiance ---------------------------------------------------

#' Get a registered irradiance calibration
#'
#' @param name Calibration id.
#'
#' @return Calibration tibble.
#' @export
get_zeus_irrad_calibration <- function(name) {
  if (!exists(name, envir = .zeus_irrad_registry, inherits = FALSE)) {
    stop("Calibration '", name, "' is not registered.", call. = FALSE)
  }
  
  get(name, envir = .zeus_irrad_registry, inherits = FALSE)
}


# Joins Irradiance/Stimulus  ----------------------------------------------

#' Join log irradiance values to stimulus metadata
#'
#' @param meta Data frame containing `wavelength` and `stim_nd`.
#' @param calib Calibration table with `wavelength`, `stim_nd`, `log_hv`.
#'
#' @return Input metadata with `log_hv`.
#' @export
nd_to_log_irradiance <- function(meta, calib) {
  required_meta <- c("wavelength", "stim_nd")
  required_calib <- c("wavelength", "stim_nd", "log_hv")
  
  if (!all(required_meta %in% names(meta))) {
    stop(
      "`meta` must contain columns: ",
      paste(required_meta, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!all(required_calib %in% names(calib))) {
    stop(
      "`calib` must contain columns: ",
      paste(required_calib, collapse = ", "),
      call. = FALSE
    )
  }
  
  meta |>
    dplyr::mutate(wavelength = as.character(.data$wavelength)) |>
    dplyr::left_join(
      calib |>
        dplyr::mutate(wavelength = as.character(.data$wavelength)),
      by = c("wavelength", "stim_nd")
    )
}


# Trace Trough-To-Peak Logic ----------------------------------------------

#' Measure one trace using ERG trough-to-peak logic
#'
#' @param trace_df Single-trace data frame with columns `time` and `value`.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param response_window_ms Numeric length-2 vector.
#' @param stimulus_onset_ms Stimulus onset in ms.
#'
#' @return One-row tibble.
#' @export
measure_trace_window <- function(trace_df,
                                 baseline_window_ms = c(300, 400),
                                 response_window_ms = c(400, 700),
                                 stimulus_onset_ms = 400) {
  trace_df <- trace_df |>
    dplyr::mutate(time_ms = zeus_time_to_ms(.data$time)) |>
    dplyr::arrange(.data$time_ms)
  
  base_df <- trace_df |>
    dplyr::filter(
      .data$time_ms >= baseline_window_ms[1],
      .data$time_ms <= baseline_window_ms[2]
    )
  
  resp_df <- trace_df |>
    dplyr::filter(
      .data$time_ms >= response_window_ms[1],
      .data$time_ms <= response_window_ms[2]
    )
  
  if (nrow(base_df) == 0L || nrow(resp_df) == 0L) {
    return(tibble::tibble(
      baseline_mv = NA_real_,
      noise_pp_mv = NA_real_,
      response_mean_mv = NA_real_,
      response_integral_mv = NA_real_,
      amp_mv = NA_real_,
      awave_mv = NA_real_,
      trough_time_ms = NA_real_,
      peak_time_ms = NA_real_,
      trough_time_poststim_ms = NA_real_,
      peak_time_poststim_ms = NA_real_
    ))
  }
  
  baseline_mv <- mean(base_df$value, na.rm = TRUE)
  noise_pp_mv <- max(base_df$value, na.rm = TRUE) - min(base_df$value, na.rm = TRUE)
  response_mean_mv <- mean(resp_df$value, na.rm = TRUE)
  response_integral_mv <- response_mean_mv - baseline_mv
  
  min_idx <- which.min(resp_df$value)
  max_idx <- which.max(resp_df$value)
  
  min_val <- resp_df$value[min_idx]
  max_val <- resp_df$value[max_idx]
  trough_time_ms <- resp_df$time_ms[min_idx]
  peak_time_ms <- resp_df$time_ms[max_idx]
  
  amp_mv <- if (isTRUE(response_integral_mv >= 0)) {
    (max_val - min_val) - noise_pp_mv
  } else {
    (min_val - max_val) + noise_pp_mv
  }
  
  tibble::tibble(
    baseline_mv = baseline_mv,
    noise_pp_mv = noise_pp_mv,
    response_mean_mv = response_mean_mv,
    response_integral_mv = response_integral_mv,
    amp_mv = amp_mv,
    awave_mv = min_val + noise_pp_mv,
    trough_time_ms = trough_time_ms,
    peak_time_ms = peak_time_ms,
    trough_time_poststim_ms = trough_time_ms - stimulus_onset_ms,
    peak_time_poststim_ms = peak_time_ms - stimulus_onset_ms
  )
}


# Extract Measurement Table -----------------------------------------------
#' Extract Irrad-wl-Amp-style measurement table
#'
#' @param x A `zeus_stimresp` object.
#' @param calib Optional calibration table with columns
#'   `wavelength`, `stim_nd`, `log_hv`.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param response_window_ms Numeric length-2 vector.
#' @param stimulus_onset_ms Numeric scalar.
#' @param same_sign Logical; if TRUE, force amplitudes to a common sign.
#'
#' @return Tibble with one row per StimResp trace.
#' @export
extract_irrad_wl_amp <- function(x,
                                 calib = NULL,
                                 baseline_window_ms = c(300, 400),
                                 response_window_ms = c(400, 700),
                                 stimulus_onset_ms = 400,
                                 same_sign = TRUE) {
  if (!inherits(x, "zeus_stimresp")) {
    stop("`x` must be a 'zeus_stimresp' object.", call. = FALSE)
  }
  
  stimresp_70 <- x$traces_70
  
  out <- stimresp_70 |>
    dplyr::group_by(
      .data$stim_index,
      .data$protocol_id,
      .data$protocol_variant,
      .data$wavelength,
      .data$stim_nd,
      .data$stim_label
    ) |>
    dplyr::group_modify(~ measure_trace_window(
      trace_df = .x,
      baseline_window_ms = baseline_window_ms,
      response_window_ms = response_window_ms,
      stimulus_onset_ms = stimulus_onset_ms
    )) |>
    dplyr::ungroup()
  
  if (!is.null(calib)) {
    out <- nd_to_log_irradiance(out, calib = calib)
  } else {
    out <- out |>
      dplyr::mutate(log_hv = NA_real_)
  }
  
  if (isTRUE(same_sign) && nrow(out) > 0L) {
    sign_mean <- mean(out$amp_mv, na.rm = TRUE)
    if (is.finite(sign_mean)) {
      out <- out |>
        dplyr::mutate(
          amp_mv = if (sign_mean >= 0) abs(.data$amp_mv) else -abs(.data$amp_mv)
        )
    }
  }
  
  out |>
    dplyr::relocate(
      .data$stim_index,
      .data$stim_label,
      .data$wavelength,
      .data$stim_nd,
      .data$log_hv
    )
}

