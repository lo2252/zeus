# validate_c0.R
# ==============================================================================
# Validate ZEUS C0 (spectral) protocol transformation against Origin export.
#
# Run this script from the package root directory:
#   Rscript validation/validate_c0.R
# Or interactively from an R session with the working directory at the
# package root (pkgload::load_all() is called automatically).
#
# REQUIRED FILES:
#   inst/extdata/26225005.abf
#   /Volumes/LOGANUSB/origin_files/26225005_origin_export_with_d_wave.xlsx
#   Or temp_file/26225005_origin_export_with_d_wave.xlsx
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load the ZEUS package
# ------------------------------------------------------------------------------

if (!requireNamespace("pkgload", quietly = TRUE)) {
  install.packages("pkgload")
}
pkgload::load_all(quiet = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(tibble)
})

# ------------------------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------------------------

abf_path    <- file.path("inst", "extdata", "26225005.abf")
xlsx_path   <- file.path("temp_file", "26225005_origin_export_with_d_wave.xlsx")
xlsx_fallback <- file.path("/Volumes", "LOGANUSB", "origin_files", "26225005_origin_export_with_d_wave.xlsx")
sheet_name  <- "StimResp"
erg_channel <- "ERG DAM80"

if (!file.exists(xlsx_path) && file.exists(xlsx_fallback)) {
  xlsx_path <- xlsx_fallback
}

stopifnot(
  "C0 ABF not found - expected at inst/extdata/26225005.abf" =
    file.exists(abf_path),
  "C0 Origin xlsx not found - expected at temp_file/... or /Volumes/LOGANUSB/origin_files/..." =
    file.exists(xlsx_path)
)

read_origin_irrad <- function(path, protocol_tbl) {
  raw_irrad <- readxl::read_excel(path, sheet = "Irrad-wl-Amp", col_names = TRUE)
  names(raw_irrad) <- make.names(names(raw_irrad), unique = TRUE)

  origin_irrad <- raw_irrad[-c(1, 2), ] |>
    dplyr::mutate(row_id = dplyr::row_number())

  numeric_cols <- c(
    "wavelength", "stimulus.wl", "stimulus.ND", "stim.irradiance",
    "response", "noise", "dwave", "peaktime", "troughtime", "dpeakime",
    "awave"
  )

  for (nm in intersect(numeric_cols, names(origin_irrad))) {
    origin_irrad[[nm]] <- suppressWarnings(as.numeric(origin_irrad[[nm]]))
  }

  protocol_tbl <- protocol_tbl |>
    dplyr::mutate(row_id = dplyr::row_number())

  origin_irrad |>
    dplyr::filter(!is.na(.data$stimulus.ND)) |>
    dplyr::left_join(
      protocol_tbl |>
        dplyr::transmute(
          row_id,
          stim_index,
          stim_label,
          protocol_id,
          wavelength_protocol = wavelength,
          stim_nd_protocol = stim_nd
        ),
      by = "row_id"
    )
}

summarise_agreement <- function(df, val_col, ref_col, label) {
  diff_col <- df[[val_col]] - df[[ref_col]]
  cat(sprintf(
    "\n  %s vs Origin:\n    n=%d  mean_diff=%.3f  mean|diff|=%.3f  RMSE=%.3f  max|diff|=%.3f  r=%.6f\n",
    label, nrow(df),
    mean(diff_col, na.rm = TRUE),
    mean(abs(diff_col), na.rm = TRUE),
    sqrt(mean(diff_col^2, na.rm = TRUE)),
    max(abs(diff_col), na.rm = TRUE),
    suppressWarnings(cor(df[[val_col]], df[[ref_col]], use = "complete.obs"))
  ))
}

# ------------------------------------------------------------------------------
# 3. Import ZEUS C0 data
# ------------------------------------------------------------------------------

cat("\n-- Importing ZEUS C0 ABF --\n")
zeus_obj <- zeus_read_abf(abf_path, protocol = "C0")
traces_280 <- zeus_obj$traces_280 |>
  dplyr::filter(!is.na(.data$stim_index))

cat(sprintf(
  "  traces_280: %d sweeps, %d rows\n",
  dplyr::n_distinct(traces_280$sweep),
  nrow(traces_280)
))

# ------------------------------------------------------------------------------
# 4. Import Origin StimResp
# ------------------------------------------------------------------------------

cat("\n-- Importing Origin C0 StimResp --\n")
raw_sheet <- readxl::read_excel(xlsx_path, sheet = sheet_name, col_names = TRUE)

time_col <- names(raw_sheet)[1]

# Trace columns match either a spectral wavelength or 650A/650B headers
trace_cols <- names(raw_sheet)[
  stringr::str_detect(
    names(raw_sheet),
    "^(650\\.00[AB]|[0-9]+\\.00)\\s+[0-9]+\\.?[0-9]*$"
  )
]

cat(sprintf("  Origin trace columns found: %d (expected 70)\n", length(trace_cols)))
if (length(trace_cols) != 70L) {
  cat("  Column names:\n")
  print(names(raw_sheet))
  stop("Expected 70 C0 trace columns in StimResp sheet.", call. = FALSE)
}

# Canonicalize "650.00A 4.0" -> "650A 4.0", "570.00 5.0" -> "570 5.0"
canon_label <- function(x) {
  x |>
    stringr::str_replace("^650\\.00A", "650A") |>
    stringr::str_replace("^650\\.00B", "650B") |>
    stringr::str_replace("^([0-9]+)\\.00", "\\1") |>
    stringr::str_squish()
}

origin_70 <- tibble::tibble(
  stim_index = seq_along(trace_cols),
  stim_col   = trace_cols,
  stim_label = canon_label(trace_cols)
)

cat("\nOrigin 70-condition protocol:\n")
print(origin_70, n = Inf)

# ------------------------------------------------------------------------------
# 5. Build Zeus 70-label protocol and compare to Origin order
# ------------------------------------------------------------------------------

cat("\n-- Comparing C0 protocol labels to Origin --\n")

zeus_70 <- protocol_table_C0() |>
  dplyr::mutate(
    stim_label = purrr::pmap_chr(
      list(.data$wavelength, .data$stim_nd, .data$stim_index, .data$protocol_id),
      make_stimresp_label
    )
  )

protocol_compare <- origin_70 |>
  dplyr::rename(origin_label = stim_label) |>
  dplyr::left_join(
    zeus_70 |> dplyr::select("stim_index", zeus_label = stim_label),
    by = "stim_index"
  ) |>
  dplyr::mutate(match = .data$origin_label == .data$zeus_label)

n_match <- sum(protocol_compare$match, na.rm = TRUE)
cat(sprintf("  Protocol label matches: %d / 70\n", n_match))

mismatches <- dplyr::filter(protocol_compare, !.data$match)
if (nrow(mismatches) > 0L) {
  cat("\n  PROTOCOL MISMATCHES:\n")
  print(mismatches, n = Inf)
} else {
  cat("  All 70 C0 protocol labels match Origin exactly.\n")
}

# ------------------------------------------------------------------------------
# 6. Compare sweep → stim_index mapping for first and last few sweeps
# ------------------------------------------------------------------------------

cat("\n-- Sweep-to-stim mapping validation (4-rep protocol) --\n")

sweep_tbl <- expand_protocol_repeats(protocol_table_C0(), repeats_per_stim = 4L) |>
  dplyr::mutate(
    sweep = .data$sweep_in_protocol,
    zeus_label = purrr::pmap_chr(
      list(.data$wavelength, .data$stim_nd, .data$stim_index, .data$protocol_id),
      make_stimresp_label
    )
  )

origin_280 <- origin_70 |>
  tidyr::uncount(weights = 4L, .id = "tech_rep") |>
  dplyr::mutate(sweep = dplyr::row_number())

sweep_compare <- origin_280 |>
  dplyr::rename(origin_label = stim_label) |>
  dplyr::left_join(
    sweep_tbl |> dplyr::select("sweep", "zeus_label"),
    by = "sweep"
  ) |>
  dplyr::mutate(match = .data$origin_label == .data$zeus_label)

n_sweep_match <- sum(sweep_compare$match, na.rm = TRUE)
cat(sprintf("  280-sweep label matches: %d / 280\n", n_sweep_match))

sweep_mismatches <- dplyr::filter(sweep_compare, !.data$match)
if (nrow(sweep_mismatches) > 0L) {
  cat("\n  SWEEP MISMATCHES:\n")
  print(sweep_mismatches, n = Inf)
} else {
  cat("  All 280 sweep labels match Origin exactly.\n")
}

# ------------------------------------------------------------------------------
# 7. Mean trace comparison: ZEUS vs Origin
#
# Origin StimResp contains BASELINE-CORRECTED, AVERAGED traces.
# The correct ZEUS counterpart is traces_70:
#   - value_raw : baseline-corrected, averaged across reps (no smoothing) — best
#                 for a direct numerical comparison to Origin StimResp.
#   - value     : same but additionally Savitzky-Golay smoothed — shown as
#                 secondary comparison to quantify smoothing effect.
# Using raw traces_280 (per-sweep, no baseline correction) would yield poor
# agreement because the DC offset is not removed.
# ------------------------------------------------------------------------------

cat("\n-- Mean trace comparison --\n")

origin_long <- raw_sheet |>
  dplyr::select(dplyr::all_of(c(time_col, trace_cols))) |>
  dplyr::rename(time_raw = dplyr::all_of(time_col)) |>
  dplyr::mutate(time_ms = suppressWarnings(as.numeric(.data$time_raw))) |>
  dplyr::filter(!is.na(.data$time_ms)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(trace_cols),
    names_to  = "stim_col",
    values_to = "origin_value"
  ) |>
  dplyr::mutate(
    origin_value = suppressWarnings(as.numeric(.data$origin_value)),
    stim_label   = canon_label(.data$stim_col)
  ) |>
  dplyr::select("stim_label", "time_ms", "origin_value")

# ZEUS processed traces (baseline-corrected, averaged, pre-smoothing)
traces_70 <- zeus_obj$traces_70

if (!"value_raw" %in% names(traces_70)) {
  # Fallback: value_raw not present (older build without the fix); use value
  traces_70 <- traces_70 |> dplyr::mutate(value_raw = .data$value)
  cat("  NOTE: value_raw not found in traces_70 - using smoothed value instead.\n")
}

# Aggregate by stim_label so that any duplicate stim_labels in traces_70
# (e.g., same ND level repeated across wavelength groups) are collapsed to a
# single row per (stim_label, time_ms) — matching what Origin exports.
# For C0 all 70 stim_labels are unique, so this is a no-op; it is included
# for consistency with validate_c1.R where it is essential.
zeus_processed <- traces_70 |>
  dplyr::select("stim_label", "time_ms",
                zeus_raw = "value_raw",
                zeus_smoothed = "value") |>
  dplyr::group_by(.data$stim_label, .data$time_ms) |>
  dplyr::summarise(
    zeus_raw      = mean(.data$zeus_raw,      na.rm = TRUE),
    zeus_smoothed = mean(.data$zeus_smoothed, na.rm = TRUE),
    .groups = "drop"
  )

wave_compare <- zeus_processed |>
  dplyr::inner_join(origin_long, by = c("stim_label", "time_ms")) |>
  dplyr::mutate(
    diff_raw      = .data$zeus_raw - .data$origin_value,
    abs_diff_raw  = abs(.data$diff_raw),
    diff_smoothed = .data$zeus_smoothed - .data$origin_value,
    abs_diff_smoothed = abs(.data$diff_smoothed)
  )

cat(sprintf("  Matched trace rows: %d\n", nrow(wave_compare)))
cat(sprintf("  Unique labels matched: %d / %d ZEUS, %d Origin\n",
            dplyr::n_distinct(wave_compare$stim_label),
            dplyr::n_distinct(zeus_processed$stim_label),
            dplyr::n_distinct(origin_long$stim_label)))

if (nrow(wave_compare) == 0L) {
  cat("  WARNING: No rows matched.\n")
  cat("  Possible causes:\n")
  cat("    - time_ms mismatch: check that Origin time column is in milliseconds\n")
  cat("    - stim_label mismatch: run protocol comparison above to verify\n")
} else {
  summarise_agreement <- function(df, val_col, ref_col, label) {
    diff_col <- df[[val_col]] - df[[ref_col]]
    cat(sprintf(
      "\n  %s vs Origin:\n    n=%d  mean_diff=%.3f  mean|diff|=%.3f  RMSE=%.3f  max|diff|=%.3f  r=%.6f\n",
      label, nrow(df),
      mean(diff_col, na.rm = TRUE),
      mean(abs(diff_col), na.rm = TRUE),
      sqrt(mean(diff_col^2, na.rm = TRUE)),
      max(abs(diff_col), na.rm = TRUE),
      suppressWarnings(cor(df[[val_col]], df[[ref_col]], use = "complete.obs"))
    ))
  }

  summarise_agreement(wave_compare, "zeus_raw",      "origin_value", "ZEUS unsmoothed (value_raw)")
  summarise_agreement(wave_compare, "zeus_smoothed", "origin_value", "ZEUS smoothed   (value)    ")

  by_label <- wave_compare |>
    dplyr::group_by(.data$stim_label) |>
    dplyr::summarise(
      n             = dplyr::n(),
      mean_abs_diff = mean(.data$abs_diff_raw, na.rm = TRUE),
      rmse          = sqrt(mean(.data$diff_raw^2, na.rm = TRUE)),
      correlation   = suppressWarnings(
        cor(.data$zeus_raw, .data$origin_value, use = "complete.obs")
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$mean_abs_diff))

  cat("\n  Per-label agreement (unsmoothed, worst first):\n")
  print(by_label, n = 10)
}

# ------------------------------------------------------------------------------
# 8. Irrad-wl-Amp comparison: per-stimulus measurements and summary statistics
# ------------------------------------------------------------------------------

cat("\n-- Comparing C0 Irrad-wl-Amp statistics to Origin --\n")

origin_irrad <- read_origin_irrad(xlsx_path, zeus_70)

zeus_irrad <- extract_irrad_wl_amp(zeus_obj)

irrad_compare <- origin_irrad |>
  dplyr::select(
    "stim_index",
    "stim_label",
    origin_response = .data$response,
    origin_awave = .data$awave,
    origin_dwave = .data$dwave,
    origin_peaktime = .data$peaktime,
    origin_troughtime = .data$troughtime,
    origin_dpeaktime = .data$dpeakime
  ) |>
  dplyr::inner_join(
    zeus_irrad |>
      dplyr::select(
        "stim_index",
        "stim_label",
        zeus_amp = .data$amp_mv,
        zeus_response_integral = .data$response_integral_mv,
        zeus_awave = .data$awave_mv,
        zeus_dwave = .data$dwave_mv,
        zeus_peaktime = .data$peak_time_poststim_ms,
        zeus_troughtime = .data$trough_time_poststim_ms,
        zeus_dpeaktime = .data$dpeak_time_poststim_ms
      ),
    by = c("stim_index", "stim_label")
  )

cat(sprintf("  Matched Irrad-wl-Amp rows: %d\n", nrow(irrad_compare)))

if (nrow(irrad_compare) > 0L) {
  summarise_agreement(irrad_compare, "zeus_amp", "origin_response", "ZEUS B-wave amplitude (amp_mv)")
  summarise_agreement(irrad_compare, "zeus_response_integral", "origin_response", "ZEUS response integral")
  summarise_agreement(irrad_compare, "zeus_awave", "origin_awave", "ZEUS A-wave amplitude")
  summarise_agreement(irrad_compare, "zeus_dwave", "origin_dwave", "ZEUS D-wave amplitude")
  summarise_agreement(irrad_compare, "zeus_peaktime", "origin_peaktime", "ZEUS B-wave peak time (post-stim)")
  summarise_agreement(irrad_compare, "zeus_troughtime", "origin_troughtime", "ZEUS A-wave trough time (post-stim)")
  summarise_agreement(irrad_compare, "zeus_dpeaktime", "origin_dpeaktime", "ZEUS D-wave peak time (post-stim)")
}

origin_irrad_for_stats <- origin_irrad |>
  dplyr::transmute(
    protocol_id = .data$protocol_id,
    stim_index = .data$stim_index,
    stim_label = .data$stim_label,
    wavelength = as.character(.data$wavelength_protocol),
    stim_nd = .data$stim_nd_protocol,
    amp_mv = .data$response,
    awave_mv = .data$awave,
    dwave_mv = .data$dwave,
    peak_time_poststim_ms = .data$peaktime,
    trough_time_poststim_ms = .data$troughtime,
    dpeak_time_poststim_ms = .data$dpeakime
  )

zeus_stats_summary <- zeus_summarize_peak_statistics(zeus_obj)
origin_stats_summary <- zeus_summarize_peak_statistics(origin_irrad_for_stats)

stats_compare <- zeus_stats_summary$combined_export |>
  dplyr::rename_with(~ paste0("zeus_", .x), c(
    "n", "mean_peak_mv", "sd_peak_mv", "sem_peak_mv", "median_peak_mv",
    "min_peak_mv", "max_peak_mv", "mean_latency_ms", "sd_latency_ms"
  )) |>
  dplyr::inner_join(
    origin_stats_summary$combined_export |>
      dplyr::rename_with(~ paste0("origin_", .x), c(
        "n", "mean_peak_mv", "sd_peak_mv", "sem_peak_mv", "median_peak_mv",
        "min_peak_mv", "max_peak_mv", "mean_latency_ms", "sd_latency_ms"
      )),
    by = c("summary_type", "protocol_id", "grouping_variable", "grouping_value", "peak_type", "wavelength", "stim_nd")
  )

cat(sprintf("  Matched summary-stat rows: %d\n", nrow(stats_compare)))

if (nrow(stats_compare) > 0L) {
  summarise_agreement(stats_compare, "zeus_mean_peak_mv", "origin_mean_peak_mv", "Summary mean peak")
  summarise_agreement(stats_compare, "zeus_sd_peak_mv", "origin_sd_peak_mv", "Summary SD peak")
  summarise_agreement(stats_compare, "zeus_mean_latency_ms", "origin_mean_latency_ms", "Summary mean latency")
}

# ------------------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------------------

cat("\n== C0 Validation Summary ==\n")
cat(sprintf("  Protocol labels match:  %d / 70\n",  n_match))
cat(sprintf("  Sweep labels match:     %d / 280\n", n_sweep_match))
if (nrow(wave_compare) > 0L) {
  diff_raw <- wave_compare$zeus_raw - wave_compare$origin_value
  cat(sprintf(
    "  Trace correlation (unsmoothed): %.6f\n",
    suppressWarnings(cor(wave_compare$zeus_raw, wave_compare$origin_value,
                         use = "complete.obs"))
  ))
  cat(sprintf("  Mean |diff| (uV):               %.4f\n", mean(abs(diff_raw), na.rm = TRUE)))
}

results_c0 <- list(
  protocol_compare = protocol_compare,
  sweep_compare    = sweep_compare,
  origin_long      = origin_long,
  traces_70        = traces_70,
  wave_compare     = wave_compare,
  origin_irrad     = origin_irrad,
  zeus_irrad       = zeus_irrad,
  irrad_compare    = irrad_compare,
  zeus_stats_summary = zeus_stats_summary,
  origin_stats_summary = origin_stats_summary,
  stats_compare    = stats_compare
)

cat("\nDone. Results stored in 'results_c0'.\n")
