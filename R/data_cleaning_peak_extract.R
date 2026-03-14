#' Find local peaks and troughs within ERG traces
#'
#' Detects local extrema within grouped long-format ERG data by comparing each
#' point to the immediately preceding and following points. A point is labeled
#' as a peak if it is greater than both neighbors, and as a trough if it is
#' less than both neighbors.
#'
#' @param df_long A long-format ERG data frame.
#' @param time_col Unquoted time column.
#' @param value_col Unquoted response-amplitude column.
#' @param group_cols Character vector of grouping columns.
#' @param time_min Optional minimum time bound in the same units as `time_col`.
#' @param time_max Optional maximum time bound in the same units as `time_col`.
#'
#' @return A tibble of detected extrema.
#' @export
erg_find_extrema <- function(df_long,
                             time_col = time_ms,
                             value_col = value,
                             group_cols = c("sweep"),
                             time_min = NULL,
                             time_max = NULL) {

  time_col <- rlang::enquo(time_col)
  value_col <- rlang::enquo(value_col)

  missing_group_cols <- setdiff(group_cols, names(df_long))
  if (length(missing_group_cols) > 0) {
    stop(
      "These grouping columns are missing from `df_long`: ",
      paste(missing_group_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- df_long |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)), !!time_col)

  if (!is.null(time_min)) {
    df <- df |>
      dplyr::filter(!!time_col >= time_min)
  }

  if (!is.null(time_max)) {
    df <- df |>
      dplyr::filter(!!time_col <= time_max)
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
    dplyr::ungroup() |>
    dplyr::mutate(
      time_ms = !!time_col,
      amplitude_uv = !!value_col
    ) |>
    dplyr::select(
      dplyr::all_of(group_cols),
      .data$time_ms,
      .data$amplitude_uv,
      .data$point_type
    )
}

#' Creates a sweep-level ERG feature table from ZEUS long-format data using a
#' biologically ordered sequence:
#' \itemize{
#'   \item a-wave trough
#'   \item b-wave peak after the a-wave trough
#'   \item d-wave trough and d-wave peak after that trough
#' }
#'
#' This function is centered on trough-to-peak measurements and is compatible
#' with the output of `read_abf()` / `abf_as_df_long()` /
#' `add_stimulus_cols_protocol()`.
#'
#' Time windows are supplied in milliseconds, but the imported ZEUS long-format
#' file stores time in seconds; the function converts time to milliseconds
#' internally.
#'
#' B-wave values are set to `NA` when any of the following conditions are met:
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
#' D-wave peak detection is also sequence-based: the d-wave peak is chosen from
#' peaks occurring after the d-wave trough.
#'
#' @param df_long A ZEUS long-format ERG data frame. Must contain at minimum
#'   `sweep`, `time`, `value`, `stim_nd`, and `stim_irradiance_log10`.
#' @param channel_filter Optional single channel name to analyze. If `NULL`,
#'   all rows are used. If your long data contains multiple channels, you should
#'   usually supply this.
#' @param a_window Numeric length-2 vector giving the a-wave interval in ms.
#' @param b_window Numeric length-2 vector giving the b-wave interval in ms.
#' @param d_window Numeric length-2 vector giving the d-wave interval in ms.
#' @param min_b_trough_to_peak_uv Minimum allowed trough-to-peak amplitude for
#'   the b-wave candidate.
#' @param min_d_trough_to_peak_uv Minimum allowed trough-to-peak amplitude for
#'   the d-wave candidate.
#'
#' @return A tibble with one row per sweep.
#'
#' @export
erg_extract_features <- function(df_long,
                                 channel_filter = NULL,
                                 a_window = c(400, 550),
                                 b_window = c(550, 700),
                                 d_window = c(700, 1000),
                                 min_b_trough_to_peak_uv = 5,
                                 min_d_trough_to_peak_uv = 5) {

  required_cols <- c(
    "sweep", "time", "value", "stim_nd", "stim_irradiance_log10"
  )

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
      stop(
        "`channel_filter` was supplied, but `df_long` has no `channel` column.",
        call. = FALSE
      )
    }

    df_work <- df_work |>
      dplyr::filter(.data$channel == channel_filter)

    if (nrow(df_work) == 0) {
      stop(
        "No rows remained after filtering to channel = '",
        channel_filter,
        "'.",
        call. = FALSE
      )
    }
  }

  df_work <- df_work |>
    dplyr::mutate(
      time_ms = .data$time * 1000
    )

  group_keys <- c("sweep", "stim_nd", "stim_irradiance_log10")

  if ("wavelength" %in% names(df_work)) {
    group_keys <- c(group_keys, "wavelength")
  }

  if ("treatment_group" %in% names(df_work)) {
    group_keys <- c(group_keys, "treatment_group")
  }

  if ("date_of_fertilization" %in% names(df_work)) {
    group_keys <- c(group_keys, "date_of_fertilization")
  }

  if ("erg_age" %in% names(df_work)) {
    group_keys <- c(group_keys, "erg_age")
  }

  # A-wave: most negative trough in a-wave interval
  a_extrema <- erg_find_extrema(
    df_long = df_work,
    time_col = time_ms,
    value_col = value,
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
      dplyr::across(dplyr::all_of(group_keys)),
      a_wave_trough_uv = .data$amplitude_uv,
      a_wave_time_ms = .data$time_ms
    )

  # B-wave: first valid positive peak after a-wave trough with minimum TTP
  b_extrema <- erg_find_extrema(
    df_long = df_work,
    time_col = time_ms,
    value_col = value,
    group_cols = group_keys,
    time_min = b_window[1],
    time_max = b_window[2]
  ) |>
    dplyr::filter(.data$point_type == "peak")

  b_candidates <- a_summary |>
    dplyr::left_join(
      b_extrema,
      by = group_keys
    ) |>
    dplyr::filter(
      .data$time_ms > .data$a_wave_time_ms,
      .data$amplitude_uv > 0
    ) |>
    dplyr::mutate(
      b_wave_trough_to_peak_uv = .data$amplitude_uv - .data$a_wave_trough_uv
    ) |>
    dplyr::filter(
      .data$b_wave_trough_to_peak_uv >= min_b_trough_to_peak_uv
    )

  b_summary <- b_candidates |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .data$time_ms, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_keys)),
      b_wave_peak_uv = .data$amplitude_uv,
      b_wave_time_ms = .data$time_ms,
      b_wave_trough_to_peak_uv = .data$b_wave_trough_to_peak_uv
    )

  # D-wave: trough first, then first valid peak after trough
  d_extrema <- erg_find_extrema(
    df_long = df_work,
    time_col = time_ms,
    value_col = value,
    group_cols = group_keys,
    time_min = d_window[1],
    time_max = d_window[2]
  )

  d_trough_summary <- d_extrema |>
    dplyr::filter(.data$point_type == "trough") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .data$amplitude_uv, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_keys)),
      d_wave_trough_uv = .data$amplitude_uv,
      d_wave_trough_time_ms = .data$time_ms
    )

  d_peak_candidates <- d_trough_summary |>
    dplyr::left_join(
      d_extrema |>
        dplyr::filter(.data$point_type == "peak"),
      by = group_keys
    ) |>
    dplyr::filter(
      .data$time_ms > .data$d_wave_trough_time_ms
    ) |>
    dplyr::mutate(
      d_wave_trough_to_peak_uv = .data$amplitude_uv - .data$d_wave_trough_uv
    ) |>
    dplyr::filter(
      .data$d_wave_trough_to_peak_uv >= min_d_trough_to_peak_uv
    )

  d_peak_summary <- d_peak_candidates |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) |>
    dplyr::slice_min(order_by = .data$time_ms, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_keys)),
      d_wave_peak_uv = .data$amplitude_uv,
      d_wave_peak_time_ms = .data$time_ms,
      d_wave_trough_to_peak_uv = .data$d_wave_trough_to_peak_uv
    )

  df_features <- a_summary |>
    dplyr::left_join(b_summary, by = group_keys) |>
    dplyr::left_join(d_trough_summary, by = group_keys) |>
    dplyr::left_join(d_peak_summary, by = group_keys) |>
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
      )
    ) |>
    dplyr::select(
      dplyr::all_of(group_keys),
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
