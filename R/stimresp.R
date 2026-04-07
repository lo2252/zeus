
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

# Zero Baseline -----------------------------------------------------------
#' Zero baseline within each trace
#'
#' @param erg_df Long ERG data.
#' @param baseline_window_ms Numeric length-2 vector.
#'
#' @return Long ERG data with baseline-zeroed `value`.
#' @export
zero_baseline_traces <- function(erg_df, baseline_window_ms = c(300, 400)) {
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


# Smoothing Traces --------------------------------------------------------

#' Smooth traces with a running mean
#'
#' @param erg_df Long ERG data.
#' @param n Integer window size.
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


# Averaging Replicates ----------------------------------------------------

#' Average technical replicates into StimResp traces
#'
#' @param erg_labeled Long ERG data with protocol labels.
#'
#' @return 70-trace StimResp tibble.
#' @export
average_technical_replicates <- function(erg_labeled) {
  erg_labeled |>
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
      .groups = "drop"
    ) |>
    dplyr::mutate(
      stim_label = purrr::pmap_chr(
        list(.data$wavelength, .data$stim_nd, .data$stim_index, .data$protocol_id),
        make_stimresp_label
      )
    )
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


# Build StimResp ----------------------------------------------------------

#' Build a StimResp object from ABF data
#'
#' This function:
#' 1. imports raw ABF data,
#' 2. maps raw sweeps to a registered protocol,
#' 3. optional baseline-zeroing and smoothing,
#' 4. optional excluding noisy technical replicates,
#' 5. averages technical replicates into 70 StimResp traces.
#'
#' @param x File path or `zeus_abf_raw` object.
#' @param protocol Protocol id. ("C0" or "C1".)
#' @param erg_channel ERG channel value in `abf_as_df_long()`.
#' @param pc_channel Photocell channel value in `abf_as_df_long()`.
#' @param repeats_per_stim Technical replicates per stimulus.
#' @param expected_stim Number of stimulus conditions.
#' @param exclude_noisy Logical.
#' @param noise_threshold Numeric SD-ratio threshold.
#' @param zero_baseline Logical.
#' @param baseline_window_ms Numeric length-2 vector.
#' @param smooth_n Running mean window size.
#'
#' @return Object of class "zeus_stimresp".
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
                           smooth_n = 1L) {
  protocol <- match.arg(protocol)

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

  if (zero_baseline) {
    erg_labeled <- zero_baseline_traces(
      erg_labeled,
      baseline_window_ms = baseline_window_ms
    )
  }

  if (smooth_n > 1L) {
    erg_labeled <- smooth_stimresp_traces(erg_labeled, n = smooth_n)
  }

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

  photocell <- pc_df |>
    dplyr::filter(.data$sweep == min(.data$sweep, na.rm = TRUE)) |>
    dplyr::mutate(time_ms = zeus_time_to_ms(.data$time))

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
      white_ir_means = white_ir_means,
      noisy_flags = noisy_flags,
      settings = list(
        protocol = protocol,
        erg_channel = erg_channel,
        pc_channel = pc_channel,
        repeats_per_stim = repeats_per_stim,
        exclude_noisy = exclude_noisy,
        noise_threshold = noise_threshold,
        zero_baseline = zero_baseline,
        baseline_window_ms = baseline_window_ms,
        smooth_n = smooth_n
      )
    ),
    class = c("zeus_stimresp", "list")
  )
}

