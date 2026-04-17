# Peak Extraction Statistics -----------------------------------------------

.zeus_peak_feature_table <- function(
    x,
    calib = NULL,
    baseline_window_ms = NULL,
    response_window_ms = NULL,
    stimulus_onset_ms = 400,
    same_sign = TRUE,
    time_reference = c("absolute", "stimulus"),
    trough_search_window_ms = NULL,
    peak_search_window_ms = NULL,
    dwave_baseline_window_ms = NULL,
    dwave_window_ms = NULL,
    dwave_trough_search_window_ms = NULL,
    dwave_peak_search_window_ms = NULL
) {
  time_reference <- match.arg(time_reference)

  if (inherits(x, "zeus_stimresp") || inherits(x, "zeus_abf")) {
    return(
      extract_irrad_wl_amp(
        x = x,
        calib = calib,
        baseline_window_ms = baseline_window_ms,
        response_window_ms = response_window_ms,
        stimulus_onset_ms = stimulus_onset_ms,
        same_sign = same_sign,
        time_reference = time_reference,
        trough_search_window_ms = trough_search_window_ms,
        peak_search_window_ms = peak_search_window_ms,
        dwave_baseline_window_ms = dwave_baseline_window_ms,
        dwave_window_ms = dwave_window_ms,
        dwave_trough_search_window_ms = dwave_trough_search_window_ms,
        dwave_peak_search_window_ms = dwave_peak_search_window_ms
      )
    )
  }

  if (!is.data.frame(x)) {
    stop(
      "`x` must be a ZEUS object or a data frame returned by `extract_irrad_wl_amp()`.",
      call. = FALSE
    )
  }

  required_cols <- c("amp_mv", "awave_mv", "dwave_mv")
  missing_cols <- setdiff(required_cols, names(x))

  if (length(missing_cols) > 0L) {
    stop(
      "`x` is missing required peak columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  x
}

.zeus_peak_long_table <- function(features_df) {
  if (!("protocol_id" %in% names(features_df))) {
    features_df <- features_df |>
      dplyr::mutate(protocol_id = NA_character_)
  }

  if (!("wavelength" %in% names(features_df))) {
    features_df <- features_df |>
      dplyr::mutate(wavelength = NA_character_)
  }

  if (!("stim_nd" %in% names(features_df))) {
    features_df <- features_df |>
      dplyr::mutate(stim_nd = NA_real_)
  }

  if (!("stim_label" %in% names(features_df))) {
    features_df <- features_df |>
      dplyr::mutate(stim_label = NA_character_)
  }

  features_df |>
    tidyr::pivot_longer(
      cols = c("awave_mv", "amp_mv", "dwave_mv"),
      names_to = "peak_metric",
      values_to = "peak_value_mv"
    ) |>
    dplyr::mutate(
      latency_ms = dplyr::case_match(
        .data$peak_metric,
        "awave_mv" ~ .data$trough_time_poststim_ms,
        "amp_mv" ~ .data$peak_time_poststim_ms,
        "dwave_mv" ~ .data$dpeak_time_poststim_ms,
        .default = NA_real_
      ),
      peak_type = dplyr::case_match(
        .data$peak_metric,
        "awave_mv" ~ "A-wave",
        "amp_mv" ~ "B-wave",
        "dwave_mv" ~ "D-wave",
        .default = .data$peak_metric
      ),
      wavelength = as.character(.data$wavelength),
      wavelength_label = dplyr::if_else(
        !is.na(.data$stim_label) & nzchar(.data$stim_label),
        stringr::str_extract(.data$stim_label, "^\\S+"),
        .data$wavelength
      )
    )
}

.zeus_peak_summary <- function(df, group_cols) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n = sum(!is.na(.data$peak_value_mv)),
      mean_peak_mv = mean(.data$peak_value_mv, na.rm = TRUE),
      sd_peak_mv = stats::sd(.data$peak_value_mv, na.rm = TRUE),
      sem_peak_mv = dplyr::if_else(n > 0, sd_peak_mv / sqrt(n), NA_real_),
      median_peak_mv = stats::median(.data$peak_value_mv, na.rm = TRUE),
      min_peak_mv = if (all(is.na(.data$peak_value_mv))) NA_real_ else min(.data$peak_value_mv, na.rm = TRUE),
      max_peak_mv = if (all(is.na(.data$peak_value_mv))) NA_real_ else max(.data$peak_value_mv, na.rm = TRUE),
      mean_latency_ms = mean(.data$latency_ms, na.rm = TRUE),
      sd_latency_ms = stats::sd(.data$latency_ms, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(
        c(
          "mean_peak_mv", "sd_peak_mv", "sem_peak_mv", "median_peak_mv",
          "min_peak_mv", "max_peak_mv", "mean_latency_ms", "sd_latency_ms"
        ),
        ~ dplyr::if_else(is.nan(.x) | is.infinite(.x), NA_real_, .x)
      )
    )
}

.zeus_peak_combined_export <- function(key_statistics, by_nd, by_wavelength) {
  key_export <- key_statistics |>
    dplyr::mutate(
      summary_type = "key_statistics",
      grouping_variable = "overall",
      grouping_value = NA_character_,
      wavelength = NA_character_,
      stim_nd = NA_real_
    )

  nd_export <- by_nd |>
    dplyr::mutate(
      summary_type = "by_nd",
      grouping_variable = "stim_nd",
      grouping_value = as.character(.data$stim_nd),
      wavelength = NA_character_
    )

  wavelength_export <- by_wavelength |>
    dplyr::mutate(
      summary_type = "by_wavelength",
      grouping_variable = "wavelength",
      grouping_value = .data$wavelength,
      stim_nd = NA_real_
    )

  dplyr::bind_rows(
    key_export,
    nd_export,
    wavelength_export
  ) |>
    dplyr::select(
      .data$summary_type,
      .data$protocol_id,
      .data$grouping_variable,
      .data$grouping_value,
      .data$peak_type,
      .data$wavelength,
      .data$stim_nd,
      .data$n,
      .data$mean_peak_mv,
      .data$sd_peak_mv,
      .data$sem_peak_mv,
      .data$median_peak_mv,
      .data$min_peak_mv,
      .data$max_peak_mv,
      .data$mean_latency_ms,
      .data$sd_latency_ms
    )
}

#' Average A-, B-, and D-wave peaks by neutral density
#'
#' Summarizes peak amplitudes extracted by [extract_irrad_wl_amp()] across
#' neutral-density (ND) levels. This is suitable for both C0 and C1 protocols.
#' The output includes the mean, standard deviation, standard error, and
#' latency summaries for each wave type.
#'
#' @inheritParams extract_irrad_wl_amp
#' @param x A ZEUS object or a data frame returned by [extract_irrad_wl_amp()].
#'
#' @return Tibble with one row per protocol × ND × wave type.
#' @export
zeus_avg_peak_by_nd <- function(
    x,
    calib = NULL,
    baseline_window_ms = NULL,
    response_window_ms = NULL,
    stimulus_onset_ms = 400,
    same_sign = TRUE,
    time_reference = c("absolute", "stimulus"),
    trough_search_window_ms = NULL,
    peak_search_window_ms = NULL,
    dwave_baseline_window_ms = NULL,
    dwave_window_ms = NULL,
    dwave_trough_search_window_ms = NULL,
    dwave_peak_search_window_ms = NULL
) {
  features_df <- .zeus_peak_feature_table(
    x = x,
    calib = calib,
    baseline_window_ms = baseline_window_ms,
    response_window_ms = response_window_ms,
    stimulus_onset_ms = stimulus_onset_ms,
    same_sign = same_sign,
    time_reference = time_reference,
    trough_search_window_ms = trough_search_window_ms,
    peak_search_window_ms = peak_search_window_ms,
    dwave_baseline_window_ms = dwave_baseline_window_ms,
    dwave_window_ms = dwave_window_ms,
    dwave_trough_search_window_ms = dwave_trough_search_window_ms,
    dwave_peak_search_window_ms = dwave_peak_search_window_ms
  )

  .zeus_peak_long_table(features_df) |>
    dplyr::filter(!is.na(.data$stim_nd)) |>
    .zeus_peak_summary(group_cols = c("protocol_id", "stim_nd", "peak_type")) |>
    dplyr::arrange(.data$protocol_id, dplyr::desc(.data$stim_nd), .data$peak_type)
}

#' Average A-, B-, and D-wave peaks by wavelength
#'
#' Summarizes peak amplitudes across wavelength blocks for spectral C0 data.
#' The output includes the mean, standard deviation, standard error, and
#' latency summaries for each wave type. White-light C1 recordings are excluded
#' from the wavelength summary because they do not vary by wavelength.
#'
#' @inheritParams extract_irrad_wl_amp
#' @param x A ZEUS object or a data frame returned by [extract_irrad_wl_amp()].
#'
#' @return Tibble with one row per protocol × wavelength × wave type.
#' @export
zeus_avg_peak_by_wavelength <- function(
    x,
    calib = NULL,
    baseline_window_ms = NULL,
    response_window_ms = NULL,
    stimulus_onset_ms = 400,
    same_sign = TRUE,
    time_reference = c("absolute", "stimulus"),
    trough_search_window_ms = NULL,
    peak_search_window_ms = NULL,
    dwave_baseline_window_ms = NULL,
    dwave_window_ms = NULL,
    dwave_trough_search_window_ms = NULL,
    dwave_peak_search_window_ms = NULL
) {
  features_df <- .zeus_peak_feature_table(
    x = x,
    calib = calib,
    baseline_window_ms = baseline_window_ms,
    response_window_ms = response_window_ms,
    stimulus_onset_ms = stimulus_onset_ms,
    same_sign = same_sign,
    time_reference = time_reference,
    trough_search_window_ms = trough_search_window_ms,
    peak_search_window_ms = peak_search_window_ms,
    dwave_baseline_window_ms = dwave_baseline_window_ms,
    dwave_window_ms = dwave_window_ms,
    dwave_trough_search_window_ms = dwave_trough_search_window_ms,
    dwave_peak_search_window_ms = dwave_peak_search_window_ms
  )

  .zeus_peak_long_table(features_df) |>
    dplyr::filter(
      !is.na(.data$wavelength_label),
      nzchar(.data$wavelength_label),
      .data$wavelength_label != "White"
    ) |>
    .zeus_peak_summary(group_cols = c("protocol_id", "wavelength_label", "peak_type")) |>
    dplyr::rename(wavelength = .data$wavelength_label) |>
    dplyr::arrange(.data$protocol_id, .data$wavelength, .data$peak_type)
}

#' Summarize key ZEUS peak statistics for export and display
#'
#' Builds a report-ready summary of A-, B-, and D-wave peak measurements using
#' [extract_irrad_wl_amp()]. The result includes an overall key-statistics
#' table, ND summaries for both C0 and C1, and wavelength summaries for C0.
#'
#' @inheritParams extract_irrad_wl_amp
#' @param x A ZEUS object or a data frame returned by [extract_irrad_wl_amp()].
#'
#' @return Named list with `key_statistics`, `by_nd`, `by_wavelength`, and
#'   `combined_export`.
#' @export
zeus_summarize_peak_statistics <- function(
    x,
    calib = NULL,
    baseline_window_ms = NULL,
    response_window_ms = NULL,
    stimulus_onset_ms = 400,
    same_sign = TRUE,
    time_reference = c("absolute", "stimulus"),
    trough_search_window_ms = NULL,
    peak_search_window_ms = NULL,
    dwave_baseline_window_ms = NULL,
    dwave_window_ms = NULL,
    dwave_trough_search_window_ms = NULL,
    dwave_peak_search_window_ms = NULL
) {
  features_df <- .zeus_peak_feature_table(
    x = x,
    calib = calib,
    baseline_window_ms = baseline_window_ms,
    response_window_ms = response_window_ms,
    stimulus_onset_ms = stimulus_onset_ms,
    same_sign = same_sign,
    time_reference = time_reference,
    trough_search_window_ms = trough_search_window_ms,
    peak_search_window_ms = peak_search_window_ms,
    dwave_baseline_window_ms = dwave_baseline_window_ms,
    dwave_window_ms = dwave_window_ms,
    dwave_trough_search_window_ms = dwave_trough_search_window_ms,
    dwave_peak_search_window_ms = dwave_peak_search_window_ms
  )

  peak_long <- .zeus_peak_long_table(features_df)

  key_statistics <- peak_long |>
    .zeus_peak_summary(group_cols = c("protocol_id", "peak_type")) |>
    dplyr::arrange(.data$protocol_id, .data$peak_type)

  by_nd <- zeus_avg_peak_by_nd(features_df)
  by_wavelength <- zeus_avg_peak_by_wavelength(features_df)
  combined_export <- .zeus_peak_combined_export(
    key_statistics = key_statistics,
    by_nd = by_nd,
    by_wavelength = by_wavelength
  )

  list(
    key_statistics = key_statistics,
    by_nd = by_nd,
    by_wavelength = by_wavelength,
    combined_export = combined_export
  )
}
