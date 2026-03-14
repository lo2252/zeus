#' Find local peaks and troughs within ERG traces
#'
#' @description
#'Detects local extrema within grouped long-format ERG data by comparing each
#' point to the immediately preceding and following points. A point is labeled
#' as a peak if it is greater than both neighbors, and as a trough if it is
#' less than both neighbors
#'
#' @param df_long A long-format data frame.
#' @param time_col Unquoted column name containing time values in milliseconds.
#' @param value_col Unquoted column name containing response amplitudes.
#' @param group_cols Character vector of column names used to group traces.
#' @param time_min Optional minimum time bound.
#' @param time_max Optional maximum time bound.
#'
#' @return A tibble of local extrema.
#' @export

erg_find_extrema <- function(df_long,
                             time_col = ms,
                             value_col = response_uv,
                             group_cols = c("sweep"),
                             time_min = NULL,
                             time_max = NULL) {
  time_col <- rlang::enquo(time_col)
  value_col <- rlang::enquo(value_col)

  df <- df_long |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)), !!time_col)

  if (!is.null(time_min)) {
    df <- dplyr::filter(df, !!time_col >= time_min)
  }
  if (!is.null(time_max)) {
    df <- dplyr::filter(df, !!time_col <= time_max)
  }

  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::mutate(
      prev_val = dplyr::lag(!!value_col),
      next_val = dplyr::lead(!!value_col),
      point_type = dplyr::case_when(
        !is.na(.data$prev_val) & !is.na(.data$next_val) &
          (!!value_col > .data$prev_val) & (!!value_col > .data$next_val) ~ "peak",
        !is.na(.data$prev_val) & !is.na(.data$next_val) &
          (!!value_col < .data$prev_val) & (!!value_col < .data$next_val) ~ "trough",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$point_type)) |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_cols)),
      time_ms = !!time_col,
      amplitude_uv = !!value_col,
      point_type = .data$point_type
    ) |>
    dplyr::ungroup()
}

#' Extract sequence-based ERG trough-to-peak features
#'
#' Creates a sweep-level ERG feature table using a biologically ordered sequence:
#' a-wave trough, then b-wave peak after the a-wave, then d-wave trough and peak
#' in the d-wave interval.
#'
#' The function is centered on trough-to-peak measurements. It returns:
#' \itemize{
#'  \item a-wave trough amplitude and time
#'   \item b-wave peak amplitude and time
#'   \item b-wave trough-to-peak amplitude (a-wave trough to b-wave peak)
#'   \item d-wave trough amplitude and time
#'   \item d-wave peak amplitude and time
#'   \item d-wave trough-to-peak amplitude
#' }
#'
#'  B-wave values are set to `NA` when any of the following conditions are met:
#' \itemize{
#'   \item b-wave peak amplitude is negative
#'   \item a-wave time occurs after b-wave time
#'   \item b-wave time is negative
#'   \item b-wave trough-to-peak amplitude is smaller than
#'         `min_b_trough_to_peak_uv`
#' }
#'
#' To reduce noisy false-positive b-waves, the b-wave is chosen as the first
#' peak after the a-wave trough that is positive and meets the minimum
#' trough-to-peak criterion.
#'
#' If `stim_nd` and `stim_irradiance` are not already present in `df_long`,
#' they may be supplied through `stim_calib` and joined by `sweep`.
#'
#' @param df_long A long-format data frame containing at minimum
#'   `sweep`, `ms`, and `response_uv`.
#' @param stim_calib Optional data frame containing `sweep`, `stim_nd`,
#'   and `stim_irradiance`.
#' @param a_window Numeric length-2 vector for the a-wave interval.
#' @param b_window Numeric length-2 vector for the b-wave interval.
#' @param d_window Numeric length-2 vector for the d-wave interval.
#' @param min_b_trough_to_peak_uv Minimum allowed trough-to-peak amplitude for
#'   the b-wave candidate.
#'
#' @return A tibble with one row per sweep.
#'
#' @examples
#' \dontrun{
#' df_features <- erg_extract_features(
#'   df_long = df_long_stim,
#'   min_b_trough_to_peak_uv = 5)
#' }
#'
#' @export
erg_extract_features <- function(df_long,
                                 stim_calib = NULL,
                                 a_window = c(400, 550),
                                 b_window = c(550, 700),
                                 d_window = c(700, 1000),
                                 min_b_trough_to_peak_uv = 5) {

  df_joined <- df_long

  if (!all(c("stim_nd", "stim_irradiance") %in% names(df_joined))) {
    if (is.null(stim_calib)) {
      stop(
        "df_long must contain 'stim_nd' and 'stim_irradiance', or ",
        "'stim_calib' must be supplied.",
        call. = FALSE
      )
    }

    df_joined <- dplyr::left_join(df_joined, stim_calib, by = "sweep")
  }

  group_key <- c("sweep", "stim_nd", "stim_irradiance")

  # A-wave: main trough in the a-wave interval
  a_extrema <- erg_find_extrema(
    df_long = df_joined,
    time_col = ms,
    value_col = response_uv,
    group_cols = group_keys,
    time_min = a_window[1],
    time_max = a_window[2]
  ) |>
    dplyr::filter(.data$point_type == "trough")

  a_summary <- a_extrema |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .data$amplitude_uv, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      .data$sweep,
      .data$stim_nd,
      .data$stim_irradiance,
      a_wave_trough_uv = .data$amplitude_uv,
      a_wave_time_ms = .data$time_ms
    )

  # B-wave: first valid positive peak after a-wave that meets the
  # minimum trough-to-peak threshold
  b_extrema <- erg_find_extrema(
    df_long = df_joined,
    time_col = ms,
    value_col = response_uv,
    group_cols = group_keys,
    time_min = b_window[1],
    time_max = b_window[2]) |>
    dplyr::filter(.data$point_type == "peak")

  b_candidates <- a_summary |>
    dplyr::left_join(
      b_extrema,
      by = c("sweep", "stim_nd", "stim_irradiance")) |>

    dplyr::filter(
      .data$time_ms > .data$a_wave_time_ms,
      .data$amplitude_uv > 0) |>

    dplyr::mutate(
      b_wave_trough_to_peak_uv = .data$amplitude_uv - .data$a_wave_trough_uv) |>

    dplyr::filter(
      .data$b_wave_trough_to_peak_uv >= min_b_trough_to_peak_uv)

  b_summary <- b_candidates |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .date$time_ms, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      .data$sweep,
      .data$stim_nd,
      .data$stim_irradiance,
      b_wave_peak_uv = .data$amplitude_uv,
      b_wave_time_ms = .data$time_ms,
      b_wave_trough_to_peak_uv
    )

  # D-wave: Trough and peak within the d-wave interval
  d_extrema <- erg_find_extrema(
    df_long = df_joined,
    time_col = ms,
    value_col = response_uv,
    group_cols = group_keys,
    time_min = d_window[1],
    time_max = d_window[2])

  d_trough_summary <- d_extrema |>
    dplyr::filter(.data$point_type == "trough") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .data$amplitude_uv, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      .data$sweep,
      .data$stim_nd,
      .data$stim_irradiance,
      d_wave_trough_uv = .data$amplitude_uv,
      d_wave_trough_time_ms = .data$time_ms
    )

  d_peak_summary <- d_extrema |>
    dplyr::filter(.data$point_type == "peak") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_max(order_by = .data$amplitude_uv, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      .data$sweep,
      .data$stim_nd,
      .data$stim_irradiance,
      d_wave_peak_uv = .data$amplitude_uv,
      d_wave_peak_time_ms = .data$time_ms
    )

  df_features <- a_summary |>
    dplyr::left_join(b_summary, by = group_keys) %>%
    dplyr::left_join(d_trough_summary, by = group_keys) %>%
    dplyr::left_join(d_peak_summary, by = group_keys) %>%
    dplyr::mutate(
      b_wave_invalid =
        is.na(.data$b_wave_peak_uv) |
        is.na(.data$b_wave_time_ms) |
        is.na(.data$a_wave_time_ms) |
        is.na(.data$b_wave_trough_to_peak_uv) |
        .data$b_wave_peak_uv < 0 |
        .data$a_wave_time_ms > .data$b_wave_time_ms |
        .data$b_wave_time_ms < 0 |
        .data$b_wave_trough_to_peak_uv < min_b_trough_to_peak_uv,

      b_wave_peak_uv = dplyr::if_else(
        .data$b_wave_invalid,
        NA_real_,
        .data$b_wave_peak_uv
      ),
      b_wave_time_ms = dplyr::if_else(
        .data$b_wave_invalid,
        NA_real_,
        .data$b_wave_time_ms
      ),
      b_wave_trough_to_peak_uv = dplyr::if_else(
        .data$b_wave_invalid,
        NA_real_,
        .data$b_wave_trough_to_peak_uv
      ),

      d_wave_trough_to_peak_uv = .data$d_wave_peak_uv - .data$d_wave_trough_uv
    ) %>%
    dplyr::select(
      .data$sweep,
      .data$stim_nd,
      .data$stim_irradiance,
      .data$a_wave_trough_uv,
      .data$a_wave_time_ms,
      .data$b_wave_peak_uv,
      .data$b_wave_time_ms,
      .data$b_wave_trough_to_peak_uv,
      .data$d_wave_trough_uv,
      .data$d_wave_trough_time_ms,
      .data$d_wave_peak_uv,
      .data$d_wave_peak_time_ms,
      .data$d_wave_trough_to_peak_uv
    )

  df_features
}


