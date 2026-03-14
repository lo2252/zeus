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




