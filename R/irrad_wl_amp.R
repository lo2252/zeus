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
#' @param trace_df Single-trace data frame.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param response_window_ms Numeric length-2 vector.
#' @param stimulus_onset_ms Numeric scalar used only for absolute-time fallback.
#' @param time_reference One of `"absolute"` or `"stimulus"`.
#' @param trough_search_window_ms Numeric length-2 vector giving the window in
#'   which to search for the trough.
#' @param peak_search_window_ms Numeric length-2 vector giving the window in
#'   which to search for the peak after the trough.
#'
#' @return One-row tibble.
#' @export
measure_trace_window <- function(trace_df,
                                 baseline_window_ms = c(300, 400),
                                 response_window_ms = c(400, 700),
                                 stimulus_onset_ms = 400,
                                 time_reference = c("absolute", "stimulus"),
                                 trough_search_window_ms = NULL,
                                 peak_search_window_ms = NULL) {
  time_reference <- match.arg(time_reference)

  trace_df <- trace_df |>
    dplyr::mutate(
      time_ms = zeus_time_to_ms(.data$time)
    ) |>
    dplyr::arrange(.data$time_ms)

  use_relative <- identical(time_reference, "stimulus") &&
    "time_rel_ms" %in% names(trace_df) &&
    any(!is.na(trace_df$time_rel_ms))

  time_col <- if (use_relative) "time_rel_ms" else "time_ms"

  # Early trough -> peak-after-trough logic
  if (is.null(trough_search_window_ms)) {
    trough_search_window_ms <- if (use_relative) c(0, 150) else c(400, 550)
  }

  if (is.null(peak_search_window_ms)) {
    peak_search_window_ms <- if (use_relative) c(0, 200) else c(400, 600)
  }

  base_df <- trace_df |>
    dplyr::filter(
      .data[[time_col]] >= baseline_window_ms[1],
      .data[[time_col]] <= baseline_window_ms[2]
    )

  resp_df <- trace_df |>
    dplyr::filter(
      .data[[time_col]] >= response_window_ms[1],
      .data[[time_col]] <= response_window_ms[2]
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


  # Find EARLY trough (a-wave)
  trough_candidates <- trace_df |>
    dplyr::filter(
      .data[[time_col]] >= trough_search_window_ms[1],
      .data[[time_col]] <= trough_search_window_ms[2]
    ) |>
    dplyr::arrange(.data[[time_col]])

  if (nrow(trough_candidates) == 0L) {
    return(tibble::tibble(
      baseline_mv = baseline_mv,
      noise_pp_mv = noise_pp_mv,
      response_mean_mv = response_mean_mv,
      response_integral_mv = response_integral_mv,
      amp_mv = NA_real_,
      awave_mv = NA_real_,
      trough_time_ms = NA_real_,
      peak_time_ms = NA_real_,
      trough_time_poststim_ms = NA_real_,
      peak_time_poststim_ms = NA_real_
    ))
  }

  trough_row <- trough_candidates |>
    dplyr::slice_min(order_by = .data$value, n = 1, with_ties = FALSE)

  min_val <- trough_row$value[[1]]
  trough_time_current <- trough_row[[time_col]][[1]]

  # Find peak AFTER trough (b-wave style)
  peak_candidates <- trace_df |>
    dplyr::filter(
      .data[[time_col]] >= max(peak_search_window_ms[1], trough_time_current),
      .data[[time_col]] <= peak_search_window_ms[2]
    ) |>
    dplyr::arrange(.data[[time_col]])

  if (nrow(peak_candidates) == 0L) {
    return(tibble::tibble(
      baseline_mv = baseline_mv,
      noise_pp_mv = noise_pp_mv,
      response_mean_mv = response_mean_mv,
      response_integral_mv = response_integral_mv,
      amp_mv = NA_real_,
      awave_mv = min_val + noise_pp_mv,
      trough_time_ms = trough_time_current,
      peak_time_ms = NA_real_,
      trough_time_poststim_ms = if (use_relative) {
        trough_time_current
      } else {
        trough_time_current - stimulus_onset_ms
      },
      peak_time_poststim_ms = NA_real_
    ))
  }

  peak_row <- peak_candidates |>
    dplyr::slice_max(order_by = .data$value, n = 1, with_ties = FALSE)

  max_val <- peak_row$value[[1]]
  peak_time_current <- peak_row[[time_col]][[1]]

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
    trough_time_ms = trough_time_current,
    peak_time_ms = peak_time_current,
    trough_time_poststim_ms = if (use_relative) {
      trough_time_current
    } else {
      trough_time_current - stimulus_onset_ms
    },
    peak_time_poststim_ms = if (use_relative) {
      peak_time_current
    } else {
      peak_time_current - stimulus_onset_ms
    }
  )
}

# Extract Measurement Table -----------------------------------------------

#' Extract an Irrad-wl-Amp-style measurement table
#'
#' @param x A `zeus_stimresp` object.
#' @param calib Optional calibration table.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param response_window_ms Numeric length-2 vector.
#' @param stimulus_onset_ms Numeric scalar used for absolute-time fallback.
#' @param same_sign Logical.
#' @param time_reference One of `"absolute"` or `"stimulus"`.
#' @param trough_search_window_ms Numeric length-2 vector giving the window in
#'   which to search for the trough.
#' @param peak_search_window_ms Numeric length-2 vector giving the window in
#'   which to search for the peak after the trough.
#'
#' @return Tibble with one row per StimResp trace.
#' @export
extract_irrad_wl_amp <- function(x,
                                 calib = NULL,
                                 baseline_window_ms = c(300, 400),
                                 response_window_ms = c(400, 700),
                                 stimulus_onset_ms = 400,
                                 same_sign = TRUE,
                                 time_reference = c("absolute", "stimulus"),
                                 trough_search_window_ms = NULL,
                                 peak_search_window_ms = NULL) {
  time_reference <- match.arg(time_reference)

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
      stimulus_onset_ms = stimulus_onset_ms,
      time_reference = time_reference,
      trough_search_window_ms = trough_search_window_ms,
      peak_search_window_ms = peak_search_window_ms
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
