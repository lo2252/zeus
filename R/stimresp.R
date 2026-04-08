# Build Stimulus Response -------------------------------------------------

#' Build a StimResp object from ABF data
#'
#' This function:
#' 1. imports raw ABF data,
#' 2. maps raw sweeps to a registered protocol,
#' 3. optionally excludes noisy technical replicates,
#' 4. averages technical replicates into 70 StimResp traces,
#' 5. post-processes the averaged StimResp traces,
#' 6. computes optional white-light means for C1.
#'
#' Build a StimResp object from ABF data
#'
#' @param x File path or `zeus_abf_raw` object.
#' @param protocol Protocol id. (`"C0"` or `"C1"`).
#' @param erg_channel ERG channel value in `abf_as_df_long()`.
#' @param pc_channel Photocell channel value in `abf_as_df_long()`.
#' @param repeats_per_stim Technical replicates per stimulus.
#' @param expected_stim Number of stimulus conditions.
#' @param exclude_noisy Logical.
#' @param noise_threshold Numeric SD-ratio threshold.
#' @param zero_baseline Logical.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param smooth_n Running mean window size.
#' @param align_to_stimulus One of `"protocol"` or `"photocell"`.
#' @param photocell_baseline_window_ms Numeric length-2 vector for photocell baseline.
#' @param photocell_threshold_frac Fraction of photocell rise used for onset detection.
#'
#' @return Object of class `"zeus_stimresp"`.
#' @export
build_stimresp <- function(x,
                           protocol = c("C0", "C1"),
                           erg_channel = "ERG DAM80",
                           pc_channel = "Photocell",
                           repeats_per_stim = 4L,
                           expected_stim = 70L,
                           exclude_noisy = FALSE,
                           noise_threshold = 1.5,
                           zero_baseline = TRUE,
                           baseline_window_ms = c(300, 400),
                           smooth_n = 1L,
                           align_to_stimulus = c("protocol", "photocell"),
                           photocell_baseline_window_ms = c(0, 300),
                           photocell_threshold_frac = 0.5) {
  protocol <- match.arg(protocol)
  align_to_stimulus <- match.arg(align_to_stimulus)

  raw_obj <- if (inherits(x, "zeus_abf_raw")) {
    x
  } else {
    read_abf_raw(x)
  }

  df_long <- abf_as_df_long(raw_obj$raw)
  validate_long_abf_df(df_long)

  erg_df <- extract_channel_trace_df(df_long, erg_channel)
  pc_df  <- extract_channel_trace_df(df_long, pc_channel)

  sweep_ids <- sort(unique(erg_df$sweep))
  expected_sweeps <- as.integer(expected_stim) * as.integer(repeats_per_stim)

  if (length(sweep_ids) < expected_sweeps) {
    stop(
      "Expected at least ", expected_sweeps,
      " ERG sweeps, found ", length(sweep_ids), ".",
      call. = FALSE
    )
  }

  if (length(sweep_ids) > expected_sweeps) {
    warning(
      "Found ", length(sweep_ids), " ERG sweeps; using first ",
      expected_sweeps, " sweeps."
    )
    sweep_ids <- sweep_ids[seq_len(expected_sweeps)]

    erg_df <- erg_df |>
      dplyr::filter(.data$sweep %in% sweep_ids)

    pc_df <- pc_df |>
      dplyr::filter(.data$sweep %in% sweep_ids)
  }

  protocol_tbl <- switch(
    protocol,
    C0 = protocol_table_C0(),
    C1 = protocol_table_C1()
  )

  protocol_sweeps <- expand_protocol_repeats(
    protocol_tbl = protocol_tbl,
    repeats_per_stim = repeats_per_stim
  ) |>
    dplyr::mutate(
      sweep = sweep_ids,
      .before = 1
    )

  erg_labeled <- erg_df |>
    dplyr::left_join(protocol_sweeps, by = "sweep")

  if (any(is.na(erg_labeled$stim_index))) {
    stop("Protocol join failed for one or more sweeps.", call. = FALSE)
  }

  # Add stimulus onset alignment
  if (identical(align_to_stimulus, "photocell")) {
    onset_df <- compute_photocell_onsets(
      pc_df = pc_df,
      baseline_window_ms = photocell_baseline_window_ms,
      threshold_frac = photocell_threshold_frac
    )
  } else {
    onset_df <- tibble::tibble(
      sweep = sweep_ids,
      stim_onset_ms = 0
    )
  }

  erg_labeled <- add_relative_time_cols(
    df = erg_labeled,
    onset_df = onset_df,
    by = "sweep",
    onset_col = "stim_onset_ms"
  )

  noisy_flags <- NULL

  if (isTRUE(exclude_noisy)) {
    noisy_flags <- compute_noisy_trace_flags(
      erg_labeled = erg_labeled,
      response_window_ms = c(300, 1000),
      noise_threshold = noise_threshold
    )

    erg_labeled <- erg_labeled |>
      dplyr::left_join(
        noisy_flags |>
          dplyr::select(.data$sweep, .data$trace_sd, .data$sd_ratio, .data$keep),
        by = "sweep"
      ) |>
      dplyr::filter(is.na(.data$keep) | .data$keep)
  }

  stimresp_70 <- average_technical_replicates(erg_labeled)

  stimresp_70 <- postprocess_stimresp_70(
    stimresp_70 = stimresp_70,
    smooth_n = smooth_n,
    zero_baseline = zero_baseline,
    baseline_window_ms = baseline_window_ms,
    time_reference = if (identical(align_to_stimulus, "photocell")) "stimulus" else "absolute"
  )

  photocell <- pc_df |>
    dplyr::filter(.data$sweep == min(.data$sweep, na.rm = TRUE)) |>
    dplyr::mutate(
      time_ms = zeus_time_to_ms(.data$time)
    )

  white_ir_means <- NULL
  if (identical(protocol, "C1")) {
    white_ir_means <- compute_white_ir_means(stimresp_70)
  }

  structure(
    list(
      raw_abf = raw_obj,
      traces_280 = erg_labeled,
      traces_70 = stimresp_70,
      protocol_70 = protocol_tbl,
      protocol_280 = protocol_sweeps,
      photocell = photocell,
      photocell_onsets = onset_df,
      white_ir_means = white_ir_means,
      noisy_flags = noisy_flags,
      settings = list(
        protocol = protocol,
        erg_channel = erg_channel,
        pc_channel = pc_channel,
        repeats_per_stim = repeats_per_stim,
        expected_stim = expected_stim,
        exclude_noisy = exclude_noisy,
        noise_threshold = noise_threshold,
        zero_baseline = zero_baseline,
        baseline_window_ms = baseline_window_ms,
        smooth_n = smooth_n,
        align_to_stimulus = align_to_stimulus,
        photocell_baseline_window_ms = photocell_baseline_window_ms,
        photocell_threshold_frac = photocell_threshold_frac
      )
    ),
    class = c("zeus_stimresp", "list")
  )
}

# Computes Noisy Traces ---------------------------------------------------

#' Compute noisy traces within technical replicate groups
#'
#' @param erg_labeled Long ERG data with `stim_index`.
#' @param response_window_ms Numeric length-2 vector.
#' @param noise_threshold Numeric threshold on SD ratio.
#'
#' @return Tibble with one row per sweep.
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
      keep = .data$sd_ratio < noise_threshold
    ) |>
    dplyr::ungroup()
}

# Zero Baseline: raw traces -----------------------------------------------

#' Zero baseline within each raw sweep trace
#'
#' @param erg_df Long ERG data with a `sweep` column.
#' @param baseline_window_ms Numeric length-2 vector.
#'
#' @return Long ERG data with baseline-zeroed `value`.
#' @export
zero_baseline_traces <- function(erg_df, baseline_window_ms = c(300, 400)) {
  if (!"sweep" %in% names(erg_df)) {
    stop("`erg_df` must contain a `sweep` column.", call. = FALSE)
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
    dplyr::select(-.data$baseline, -.data$time_ms)
}

# Smoothing: raw traces ---------------------------------------------------

#' Smooth raw traces with a running mean
#'
#' @param erg_df Long ERG data with a `sweep` column.
#' @param n Integer window size.
#'
#' @return Long ERG data with smoothed `value`.
#' @export
smooth_stimresp_traces <- function(erg_df, n = 1L) {
  n <- as.integer(n)

  if (!"sweep" %in% names(erg_df)) {
    stop("`erg_df` must contain a `sweep` column.", call. = FALSE)
  }

  if (n <= 1L) {
    return(erg_df)
  }

  erg_df |>
    dplyr::group_by(.data$sweep) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::mutate(
      value = run_mean(.data$value, n = n)
    ) |>
    dplyr::ungroup()
}

# Averaging Replicates ----------------------------------------------------

#' Average technical replicates into StimResp traces
#'
#' @param erg_labeled Long ERG data with protocol labels.
#'
#' @return 70-trace StimResp tibble.
#' @export
average_technical_replicates <- function(erg_labeled) {
  has_onset <- "stim_onset_ms" %in% names(erg_labeled)

  out <- erg_labeled |>
    dplyr::group_by(
      .data$stim_index,
      .data$protocol_id,
      .data$protocol_variant,
      .data$stim_type,
      .data$block_index,
      .data$within_block_index,
      .data$wavelength,
      .data$stim_nd,
      .data$time
    ) |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      n_repeats_used = dplyr::n_distinct(.data$sweep),
      stim_onset_ms = if (has_onset) {
        mean(.data$stim_onset_ms, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      stim_label = purrr::pmap_chr(
        list(.data$wavelength, .data$stim_nd, .data$stim_index, .data$protocol_id),
        make_stimresp_label
      ),
      time_ms = zeus_time_to_ms(.data$time),
      time_rel_ms = .data$time_ms - .data$stim_onset_ms
    )

  out
}

# Post-process averaged StimResp traces -----------------------------------

#' Post-process averaged StimResp traces
#'
#' @param stimresp_70 A 70-trace StimResp tibble.
#' @param smooth_n Integer running-mean window size.
#' @param zero_baseline Logical.
#' @param baseline_window_ms Numeric length-2 vector in milliseconds.
#' @param time_reference One of `"absolute"` or `"stimulus"`.
#'
#' @return A post-processed StimResp tibble.
#' @export
postprocess_stimresp_70 <- function(stimresp_70,
                                    smooth_n = 1L,
                                    zero_baseline = TRUE,
                                    baseline_window_ms = c(300, 400),
                                    time_reference = c("absolute", "stimulus")) {
  time_reference <- match.arg(time_reference)

  needed <- c("stim_index", "time", "value")
  missing_cols <- setdiff(needed, names(stimresp_70))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in `stimresp_70`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- stimresp_70 |>
    dplyr::mutate(
      time_ms = zeus_time_to_ms(.data$time)
    )

  if ("stim_onset_ms" %in% names(out)) {
    out <- out |>
      dplyr::mutate(
        time_rel_ms = .data$time_ms - .data$stim_onset_ms
      )
  } else {
    out <- out |>
      dplyr::mutate(
        time_rel_ms = NA_real_
      )
  }

  smooth_n <- as.integer(smooth_n)

  if (smooth_n > 1L) {
    out <- out |>
      dplyr::group_by(.data$stim_index) |>
      dplyr::arrange(.data$time, .by_group = TRUE) |>
      dplyr::mutate(
        value = run_mean(.data$value, n = smooth_n)
      ) |>
      dplyr::ungroup()
  }

  if (isTRUE(zero_baseline)) {
    use_relative <- identical(time_reference, "stimulus") &&
      "time_rel_ms" %in% names(out) &&
      any(!is.na(out$time_rel_ms))

    if (use_relative) {
      out <- out |>
        dplyr::group_by(.data$stim_index) |>
        dplyr::mutate(
          baseline = mean(
            .data$value[
              .data$time_rel_ms >= baseline_window_ms[1] &
                .data$time_rel_ms <= baseline_window_ms[2]
            ],
            na.rm = TRUE
          ),
          value = .data$value - .data$baseline
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.data$baseline)
    } else {
      out <- out |>
        dplyr::group_by(.data$stim_index) |>
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
  }

  out
}

# C1 Mean Trace -----------------------------------------------------------

#' Compute C1 mean traces by ND across 10 runs
#'
#' @param stimresp_70 A 70-trace StimResp tibble.
#'
#' @return Tibble with 7 mean white-light intensity traces.
#' @export
compute_white_ir_means <- function(stimresp_70) {
  stimresp_70 |>
    dplyr::filter(.data$protocol_id == "C1") |>
    dplyr::group_by(.data$stim_nd, .data$time) |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      n_conditions = dplyr::n_distinct(.data$stim_index),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      wavelength = "White",
      stim_label = paste0("White ", formatC(.data$stim_nd, format = "f", digits = 1))
    )
}



# Photocell onset detection ------------------------------------------------

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

  baseline_val <- median(signal[baseline_idx], na.rm = TRUE)
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

# Compute Photocell -------------------------------------------------------

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
