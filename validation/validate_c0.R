# validate_c0.R
# ==============================================================================
# Validate ZEUS C0 (spectral) protocol transformation against Origin export.
#
# Run this script from the package root directory:
#   Rscript validation/validate_c0.R
# Or interactively from an R session with the working directory at the
# package root (pkgload::load_all() is called automatically).
#
# REQUIRED FILES (relative to package root):
#   temp_file/26225005.abf                           -- C0 ABF recording
#   temp_file/26225005_origin_export_with_d_wave.xlsx -- Origin StimResp export
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

abf_path    <- file.path("temp_file", "26225005.abf")
xlsx_path   <- file.path("temp_file", "26225005_origin_export_with_d_wave.xlsx")
sheet_name  <- "StimResp"
erg_channel <- "ERG DAM80"

stopifnot(
  "C0 ABF not found - expected at temp_file/26225005.abf" =
    file.exists(abf_path),
  "C0 Origin xlsx not found - expected at temp_file/26225005_origin_export_with_d_wave.xlsx" =
    file.exists(xlsx_path)
)

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
    zeus_70 |> dplyr::select(.data$stim_index, zeus_label = stim_label),
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
    sweep_tbl |> dplyr::select(.data$sweep, zeus_label),
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
  dplyr::select(.data$stim_label, .data$time_ms, .data$origin_value)

# ZEUS mean trace per stim_label
zeus_mean <- traces_280 |>
  dplyr::filter(.data$channel == erg_channel) |>
  dplyr::group_by(.data$stim_label, .data$time_ms) |>
  dplyr::summarise(zeus_value = mean(.data$value, na.rm = TRUE), .groups = "drop")

wave_compare <- zeus_mean |>
  dplyr::inner_join(origin_long, by = c("stim_label", "time_ms")) |>
  dplyr::mutate(
    diff      = .data$zeus_value - .data$origin_value,
    abs_diff  = abs(.data$diff),
    sq_diff   = .data$diff^2
  )

cat(sprintf("  Matched trace rows: %d\n", nrow(wave_compare)))

if (nrow(wave_compare) == 0L) {
  cat("  WARNING: No rows matched. Check time alignment or label format.\n")
} else {
  overall <- wave_compare |>
    dplyr::summarise(
      n_points     = dplyr::n(),
      mean_diff    = mean(.data$diff,     na.rm = TRUE),
      mean_abs_diff = mean(.data$abs_diff, na.rm = TRUE),
      rmse         = sqrt(mean(.data$sq_diff, na.rm = TRUE)),
      max_abs_diff = max(.data$abs_diff,  na.rm = TRUE),
      correlation  = suppressWarnings(
        cor(.data$zeus_value, .data$origin_value, use = "complete.obs")
      )
    )

  cat("\n  Overall trace agreement:\n")
  print(overall)

  by_label <- wave_compare |>
    dplyr::group_by(.data$stim_label) |>
    dplyr::summarise(
      n            = dplyr::n(),
      mean_abs_diff = mean(.data$abs_diff, na.rm = TRUE),
      rmse         = sqrt(mean(.data$sq_diff, na.rm = TRUE)),
      correlation  = suppressWarnings(
        cor(.data$zeus_value, .data$origin_value, use = "complete.obs")
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$mean_abs_diff))

  cat("\n  Per-label trace agreement (worst first):\n")
  print(by_label, n = 10)
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------

cat("\n== C0 Validation Summary ==\n")
cat(sprintf("  Protocol labels match:  %d / 70\n",  n_match))
cat(sprintf("  Sweep labels match:     %d / 280\n", n_sweep_match))
if (nrow(wave_compare) > 0L) {
  cat(sprintf(
    "  Trace correlation:      %.6f\n",
    suppressWarnings(cor(wave_compare$zeus_value, wave_compare$origin_value,
                         use = "complete.obs"))
  ))
  cat(sprintf(
    "  Mean |diff| (uV):       %.4f\n",
    mean(wave_compare$abs_diff, na.rm = TRUE)
  ))
}

results_c0 <- list(
  protocol_compare = protocol_compare,
  sweep_compare    = sweep_compare,
  wave_compare     = wave_compare
)

cat("\nDone. Results stored in 'results_c0'.\n")
