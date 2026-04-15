# Build Stimulus Response -------------------------------------------------

#' Build stimulus-response traces
#'
#' Applies the standard StimResp pipeline to labeled replicate-level ERG traces:
#' optional per-sweep baseline correction, optional noisy-replicate exclusion,
#' replicate averaging, and optional post-averaging smoothing.
#'
#' @param traces_280 A long replicate-level trace table containing at least
#'   `sweep`, `stim_index`, `stim_label`, `time_ms`, and `value`.
#' @param stimresp_zero_baseline Logical; whether to baseline-correct each sweep
#'   before averaging. Default is `TRUE`.
#' @param stimresp_baseline_window_ms Numeric length-2 vector giving the
#'   baseline window in milliseconds. Default is `c(300, 400)`.
#' @param stimresp_exclude_noisy Logical; whether to exclude noisy technical
#'   replicates before averaging. Default is `TRUE`.
#' @param stimresp_noise_window_ms Numeric length-2 vector giving the time
#'   window in milliseconds used to estimate replicate-level standard deviation.
#'   Default is `c(300, 1000)`.
#' @param stimresp_noise_ratio_cutoff Numeric cutoff applied to the within-group
#'   SD ratio. Replicates with `sd_ratio >= stimresp_noise_ratio_cutoff` are
#'   removed before averaging. Default is `1.5`.
#' @param stimresp_runmean_k Integer running-average width applied after
#'   replicate averaging. Use `1L` for no running average. Default is `16L`.
#' @param stimresp_sg_smooth Logical; whether to apply Savitzky-Golay smoothing
#'   after replicate averaging. Default is `TRUE`.
#' @param stimresp_sg_n Odd integer Savitzky-Golay window size. If even, it is
#'   incremented to the next odd integer. Default is `101L`.
#' @param stimresp_sg_p Integer Savitzky-Golay polynomial order. Default is `2L`.
#'
#' @return A list with:
#' \describe{
#'   \item{traces_70}{Averaged and optionally smoothed stimulus-response traces.}
#'   \item{stimresp_qc}{Per-sweep noisy-replicate QC summary.}
#' }
#' @keywords internal
build_stimresp <- function(traces_280,
                           stimresp_zero_baseline = TRUE,
                           stimresp_baseline_window_ms = c(300, 400),
                           stimresp_exclude_noisy = TRUE,
                           stimresp_noise_window_ms = c(300, 1000),
                           stimresp_noise_ratio_cutoff = 1.5,
                           stimresp_runmean_k = 16L,
                           stimresp_sg_smooth = TRUE,
                           stimresp_sg_n = 101L,
                           stimresp_sg_p = 2L) {
  needed <- c("sweep", "stim_index", "stim_label", "time_ms", "value")
  missing_cols <- setdiff(needed, names(traces_280))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in `traces_280`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  traces_proc <- traces_280

  if (isTRUE(stimresp_zero_baseline)) {
    traces_proc <- traces_proc |>
      dplyr::group_by(.data$sweep) |>
      dplyr::group_modify(~ baseline_one_trace(
        .x,
        baseline_window_ms = stimresp_baseline_window_ms
      )) |>
      dplyr::ungroup()
  }

  avg_obj <- average_stimresp_reps(
    df = traces_proc,
    exclude_noisy = stimresp_exclude_noisy,
    noise_window_ms = stimresp_noise_window_ms,
    noise_ratio_cutoff = stimresp_noise_ratio_cutoff
  )

  # Preserve the unsmoothed average as value_raw before any post-processing
  # so that compare_raw = TRUE in zeus_plot_mean_waveform has a value to show.
  traces_70 <- avg_obj$avg |>
    dplyr::mutate(value_raw = .data$value)

  if (isTRUE(stimresp_sg_smooth) || as.integer(stimresp_runmean_k) > 1L) {
    traces_70 <- traces_70 |>
      dplyr::group_by(.data$stim_index) |>
      dplyr::group_modify(~ smooth_one_trace(
        .x,
        runmean_k = stimresp_runmean_k,
        sg_n = if (isTRUE(stimresp_sg_smooth)) stimresp_sg_n else 1L,
        sg_p = stimresp_sg_p
      )) |>
      dplyr::ungroup()
  }

  list(
    traces_70 = traces_70,
    stimresp_qc = avg_obj$qc
  )
}

# Noisy replicate detection -----------------------------------------------

#' Compute noisy replicate flags within stimulus groups
#'
#' For each stimulus group, the replicate-level standard deviation is computed
#' within a specified response window. Each replicate SD is divided by the
#' minimum SD within the same stimulus group to form an SD ratio. Replicates
#' with ratios greater than or equal to the supplied threshold are flagged for
#' exclusion.
#'
#' @param erg_labeled Long ERG data with at least `sweep`, `time`, `value`, and
#'   `stim_index`.
#' @param response_window_ms Numeric length-2 vector giving the response window
#'   in milliseconds. Default is `c(300, 1000)`.
#' @param noise_threshold Numeric SD-ratio threshold. Replicates with
#'   `sd_ratio >= noise_threshold` are flagged as noisy. Default is `1.5`.
#'
#' @return A tibble with one row per sweep containing replicate-level noise
#'   metrics and a logical `keep` column.
#' @export
compute_noisy_trace_flags <- function(erg_labeled,
                                      response_window_ms = c(300, 1000),
                                      noise_threshold = 1.5) {
  needed <- c("sweep", "time", "value", "stim_index")
  missing_cols <- setdiff(needed, names(erg_labeled))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  erg_labeled |>
    dplyr::mutate(time_ms = zeus_time_to_ms(.data$time)) |>
    dplyr::filter(
      .data$time_ms >= response_window_ms[1],
      .data$time_ms <= response_window_ms[2]
    ) |>
    dplyr::group_by(.data$sweep, .data$stim_index) |>
    dplyr::summarise(
      trace_sd = stats::sd(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$stim_index) |>
    dplyr::mutate(
      min_sd = min(.data$trace_sd, na.rm = TRUE),
      sd_ratio = .data$trace_sd / .data$min_sd,
      keep = is.finite(.data$sd_ratio) & (.data$sd_ratio < noise_threshold)
    ) |>
    dplyr::ungroup()
}

# Zero baseline: raw traces -----------------------------------------------

#' Zero baseline within each raw sweep trace
#'
#' @param erg_df Long ERG data with `sweep`, `time`, and `value`.
#' @param baseline_window_ms Numeric length-2 vector giving the baseline window
#'   in milliseconds. Default is `c(300, 400)`.
#'
#' @return Long ERG data with baseline-zeroed `value`.
#' @export
zero_baseline_traces <- function(erg_df, baseline_window_ms = c(300, 400)) {
  needed <- c("sweep", "time", "value")
  missing_cols <- setdiff(needed, names(erg_df))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in `erg_df`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  erg_df |>
    dplyr::mutate(time_ms = zeus_time_to_ms(.data$time)) |>
    dplyr::group_by(.data$sweep) |>
    dplyr::mutate(
      baseline = mean(
        .data$value[
          .data$time_ms >= baseline_window_ms[1] &
            .data$time_ms <= baseline_window_ms[2]
        ],
        na.rm = TRUE
      ),
      value = .data$value - .data$baseline
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$baseline)
}

# Smoothing: raw traces ---------------------------------------------------

#' Smooth raw traces with a running mean
#'
#' @param erg_df Long ERG data with a `sweep` column.
#' @param n Integer running-mean window size. Use `1L` for no smoothing.
#'
#' @return Long ERG data with smoothed `value`.
#' @export
smooth_stimresp_traces <- function(erg_df, n = 1L) {
  n <- as.integer(n)

  erg_df |>
    dplyr::group_by(.data$sweep) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::mutate(value = run_mean(.data$value, n = n)) |>
    dplyr::ungroup()
}

# Photocell onset detection -----------------------------------------------

#' Detect stimulus onset from a single photocell trace
#'
#' @param time Numeric vector of time values.
#' @param signal Numeric vector of photocell values.
#' @param baseline_window_ms Numeric length-2 vector used to estimate baseline.
#' @param threshold_frac Fraction of the rise from baseline to peak used as the
#'   onset threshold.
#'
#' @return Numeric scalar giving stimulus onset in milliseconds.
#' @keywords internal
detect_photocell_onset_single <- function(time,
                                          signal,
                                          baseline_window_ms = c(0, 300),
                                          threshold_frac = 0.5) {
  time_ms <- zeus_time_to_ms(time)

  baseline_idx <- time_ms >= baseline_window_ms[1] &
    time_ms <= baseline_window_ms[2]

  if (!any(baseline_idx, na.rm = TRUE)) {
    return(NA_real_)
  }

  baseline_val <- stats::median(signal[baseline_idx], na.rm = TRUE)
  peak_val <- max(signal, na.rm = TRUE)

  if (!is.finite(baseline_val) || !is.finite(peak_val) || peak_val <= baseline_val) {
    return(NA_real_)
  }

  threshold_val <- baseline_val + threshold_frac * (peak_val - baseline_val)

  onset_idx <- which(signal >= threshold_val)[1]

  if (is.na(onset_idx) || length(onset_idx) == 0L) {
    return(NA_real_)
  }

  time_ms[onset_idx]
}

# Compute photocell onsets ------------------------------------------------

#' Compute photocell-based stimulus onset for each sweep
#'
#' @param pc_df Long photocell data containing `sweep`, `time`, and `value`.
#' @param baseline_window_ms Numeric length-2 vector used to estimate baseline.
#' @param threshold_frac Fraction of the rise from baseline to peak used as the
#'   onset threshold.
#'
#' @return Tibble with columns `sweep` and `stim_onset_ms`.
#' @export
compute_photocell_onsets <- function(pc_df,
                                     baseline_window_ms = c(0, 300),
                                     threshold_frac = 0.5) {
  needed <- c("sweep", "time", "value")
  missing_cols <- setdiff(needed, names(pc_df))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in `pc_df`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  pc_df |>
    dplyr::group_by(.data$sweep) |>
    dplyr::summarise(
      stim_onset_ms = detect_photocell_onset_single(
        time = .data$time,
        signal = .data$value,
        baseline_window_ms = baseline_window_ms,
        threshold_frac = threshold_frac
      ),
      .groups = "drop"
    )
}

# Relative time helper ----------------------------------------------------

#' Add time columns relative to stimulus onset
#'
#' @param df Data frame containing `time`.
#' @param onset_df Data frame containing onset information by grouping variable.
#' @param by Column name used to join onset values, usually `"sweep"` or
#'   `"stim_index"`.
#' @param onset_col Name of onset column. Default is `"stim_onset_ms"`.
#'
#' @return Input data with `time_ms` and `time_rel_ms`.
#' @keywords internal
add_relative_time_cols <- function(df,
                                   onset_df,
                                   by = "sweep",
                                   onset_col = "stim_onset_ms") {
  if (!"time" %in% names(df)) {
    stop("`df` must contain a `time` column.", call. = FALSE)
  }

  if (!by %in% names(df)) {
    stop("`df` must contain column `", by, "`.", call. = FALSE)
  }

  if (!by %in% names(onset_df)) {
    stop("`onset_df` must contain column `", by, "`.", call. = FALSE)
  }

  if (!onset_col %in% names(onset_df)) {
    stop("`onset_df` must contain column `", onset_col, "`.", call. = FALSE)
  }

  df |>
    dplyr::left_join(onset_df, by = by) |>
    dplyr::mutate(
      time_ms = zeus_time_to_ms(.data$time),
      time_rel_ms = .data$time_ms - .data[[onset_col]]
    )
}

# Internal helpers --------------------------------------------------------

#' Baseline-correct one trace
#'
#' @param df Data frame with columns `time_ms` and `value`.
#' @param baseline_window_ms Numeric length-2 vector giving the baseline window
#'   in milliseconds.
#'
#' @return A data frame with baseline-corrected `value`.
#' @keywords internal
baseline_one_trace <- function(df, baseline_window_ms = c(300, 400)) {
  base_vals <- df$value[
    df$time_ms >= baseline_window_ms[1] &
      df$time_ms <= baseline_window_ms[2]
  ]

  base_val <- if (length(base_vals) == 0L || all(is.na(base_vals))) {
    0
  } else {
    mean(base_vals, na.rm = TRUE)
  }

  if (!is.finite(base_val)) {
    base_val <- 0
  }

  df |>
    dplyr::mutate(value = value - base_val)
}

#' Compute replicate-level noise metrics for StimResp
#'
#' For each stimulus group, compute the standard deviation over the supplied
#' noise window for each technical replicate, normalize each replicate SD by the
#' minimum SD in that group, and flag replicates with ratios below the cutoff.
#'
#' @param df Data frame with columns `stim_index`, `sweep`, `time_ms`, and
#'   `value`.
#' @param noise_window_ms Numeric length-2 vector giving the SD window in
#'   milliseconds.
#' @param noise_ratio_cutoff Numeric cutoff. Replicates with
#'   `sd_ratio >= noise_ratio_cutoff` are dropped.
#'
#' @return A data frame with per-replicate noise metrics and a logical `keep`
#'   column.
#' @keywords internal
compute_stimresp_noise_flags <- function(df,
                                         noise_window_ms = c(300, 1000),
                                         noise_ratio_cutoff = 1.5) {
  df |>
    dplyr::group_by(.data$stim_index, .data$sweep) |>
    dplyr::summarise(
      sweep_sd = stats::sd(
        .data$value[
          .data$time_ms >= noise_window_ms[1] &
            .data$time_ms <= noise_window_ms[2]
        ],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$stim_index) |>
    dplyr::mutate(
      min_sd_in_group = min(.data$sweep_sd, na.rm = TRUE),
      sd_ratio = .data$sweep_sd / .data$min_sd_in_group,
      keep = is.finite(.data$sd_ratio) & (.data$sd_ratio < noise_ratio_cutoff)
    ) |>
    dplyr::ungroup()
}

#' Average StimResp replicates with optional noisy-replicate exclusion
#'
#' @param df Data frame with columns `stim_index`, `stim_label`, `sweep`,
#'   `time_ms`, and `value`.
#' @param exclude_noisy Logical; whether to apply noisy-replicate exclusion.
#' @param noise_window_ms Numeric length-2 vector giving the SD window in
#'   milliseconds.
#' @param noise_ratio_cutoff Numeric cutoff. Replicates with
#'   `sd_ratio >= noise_ratio_cutoff` are dropped.
#'
#' @return A list with:
#' \describe{
#'   \item{avg}{Averaged traces.}
#'   \item{qc}{Per-replicate QC summary.}
#' }
#' @keywords internal
average_stimresp_reps <- function(df,
                                  exclude_noisy = TRUE,
                                  noise_window_ms = c(300, 1000),
                                  noise_ratio_cutoff = 1.5) {
  qc <- compute_stimresp_noise_flags(
    df = df,
    noise_window_ms = noise_window_ms,
    noise_ratio_cutoff = noise_ratio_cutoff
  )

  df_use <- if (isTRUE(exclude_noisy)) {
    df |>
      dplyr::left_join(
        qc |>
          dplyr::select(.data$stim_index, .data$sweep, .data$keep),
        by = c("stim_index", "sweep")
      ) |>
      dplyr::filter(.data$keep) |>
      dplyr::select(-.data$keep)
  } else {
    df
  }

  # Preserve protocol metadata columns that are constant within each
  # stim_index (e.g., wavelength, stim_nd). Including them in group_by
  # does not change the averaging logic but retains them in the output,
  # which is critical for spectral protocols like C0.
  meta_cols <- intersect(
    c("protocol_id", "protocol_variant", "stim_type",
      "block_index", "within_block_index", "wavelength", "stim_nd"),
    names(df_use)
  )

  # Also preserve the original time column if present (time_ms is always
  # grouped on; time in seconds has a 1:1 mapping and is needed by
  # downstream post-processing and plotting).
  time_cols <- intersect(c("time", "time_ms"), names(df_use))

  group_cols <- unique(c("stim_index", "stim_label", meta_cols, time_cols))

  avg <- df_use |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      n_reps_used = dplyr::n_distinct(.data$sweep),
      .groups = "drop"
    )

  list(avg = avg, qc = qc)
}

#' Smooth one averaged trace
#'
#' @param df Data frame with column `value`.
#' @param runmean_k Integer running-mean window. Use `1L` for none.
#' @param sg_n Odd integer Savitzky-Golay window size.
#' @param sg_p Integer Savitzky-Golay polynomial order.
#'
#' @return Smoothed data frame.
#' @keywords internal
smooth_one_trace <- function(df,
                             runmean_k = 16L,
                             sg_n = 101L,
                             sg_p = 2L) {
  runmean_k <- as.integer(runmean_k)
  sg_n <- as.integer(sg_n)
  sg_p <- as.integer(sg_p)

  if (sg_n %% 2L == 0L) {
    sg_n <- sg_n + 1L
  }

  out <- df

  if (runmean_k > 1L) {
    out <- out |>
      dplyr::mutate(
        value = zoo::rollmean(.data$value, k = runmean_k, fill = "extend")
      )
  }

  if (sum(is.finite(out$value)) >= sg_n) {
    out <- out |>
      dplyr::mutate(
        value = signal::sgolayfilt(.data$value, p = sg_p, n = sg_n)
      )
  }

  out
}

# StimResp QC summary -----------------------------------------------------

#' Summarize StimResp noisy-replicate exclusion
#'
#' @param x A ZEUS object containing `stimresp_qc`.
#'
#' @return A list with:
#' \describe{
#'   \item{overall}{Overall counts of retained and excluded replicates.}
#'   \item{by_stim}{Retained and excluded replicate counts by stimulus.}
#'   \item{dropped}{A table of excluded replicates.}
#' }
#' @export
summarize_stimresp_qc <- function(x) {
  if (is.null(x$stimresp_qc)) {
    stop("`x` does not contain `stimresp_qc`.", call. = FALSE)
  }

  overall <- x$stimresp_qc |>
    dplyr::summarise(
      total_sweeps = dplyr::n(),
      kept_sweeps = sum(.data$keep, na.rm = TRUE),
      dropped_sweeps = sum(!.data$keep, na.rm = TRUE),
      drop_rate = mean(!.data$keep, na.rm = TRUE)
    )

  by_stim <- x$stimresp_qc |>
    dplyr::group_by(.data$stim_index) |>
    dplyr::summarise(
      kept = sum(.data$keep, na.rm = TRUE),
      dropped = sum(!.data$keep, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$dropped), .data$stim_index)

  dropped <- x$stimresp_qc |>
    dplyr::filter(!.data$keep) |>
    dplyr::arrange(.data$stim_index, .data$sweep)

  list(
    overall = overall,
    by_stim = by_stim,
    dropped = dropped
  )
}