# ==============================================================================
# COMPARE UPDATED C0 PROTOCOL (protocol_table_C0) TO ORIGIN FILE
#
# ASSUMPTIONS
# - Your updated protocol object name is preserved as: protocol_table_C0
# - You want to import BOTH:
#     1) ZEUS long data from your current R object / import pipeline
#     2) Origin workbook sheets
#
# WHAT THIS DOES
# 1) imports Origin StimResp
# 2) imports ZEUS long data (you plug in your ZEUS import step)
# 3) standardizes the updated protocol_table_C0
# 4) compares protocol order/labels against Origin StimResp column order
# 5) rebuilds ZEUS mean traces using the UPDATED C0 protocol
# 6) compares ZEUS traces to Origin StimResp traces
#
# REQUIRED PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(tibble)
library(zoo)

if (!exists("zeus_read_abf", mode = "function") ||
    !exists("protocol_table_C0", mode = "function")) {
  source_files <- file.path(
    "R",
    c(
      "utils.R",
      "protocols.R",
      "stimresp.R",
      "ERG_Import.R",
      "irrad_wl_amp.R",
      "extdata.R",
      "plotting_function.R"
    )
  )
  for (source_file in source_files[file.exists(source_files)]) {
    sys.source(source_file, envir = globalenv())
  }
}

# ==============================================================================
# USER INPUTS
# ==============================================================================

worksheet_path <- "../test_data/Origin_Files_with_dwave/26225005_origin_export_with_d_wave.xlsx"
stimresp_sheet <- "StimResp"
channel_keep   <- "ERG DAM80"

# ------------------------------------------------------------------------------
# ZEUS IMPORT STEP
# Replace this section with your actual ZEUS import code if needed.
#
# The final object MUST be named: df_long
# and contain at least these columns:
#   wavelength, stim_nd, time_ms, value, sweep, channel
# ------------------------------------------------------------------------------

# Example placeholder:
zeus_obj <- zeus_read_abf(
  "../test_data/Origin_Files_with_dwave/26225005.abf",
  protocol = "C0"
)
df_long <- zeus_obj$traces_280

# If df_long already exists in memory, this check will pass.
stopifnot(exists("df_long"))
stopifnot(is.data.frame(df_long))
stopifnot(exists("protocol_table_C0", mode = "function"))

# ==============================================================================
# HELPERS
# ==============================================================================

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_chr_col <- function(data, nm) {
  if (nm %in% names(data)) as.character(data[[nm]]) else rep(NA_character_, nrow(data))
}

safe_num_col <- function(data, nm) {
  if (nm %in% names(data)) suppressWarnings(as.numeric(data[[nm]])) else rep(NA_real_, nrow(data))
}

canonicalize_stim_label <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace("^650\\.00A", "650A") |>
    stringr::str_replace("^650\\.00B", "650B") |>
    stringr::str_replace("^([0-9]+)\\.00", "\\1") |>
    stringr::str_squish()
}

make_stim_label <- function(wavelength, stim_nd, stim_index = NA_integer_, protocol_id = "C0") {
  if (exists("make_stimresp_label", mode = "function")) {
    return(make_stimresp_label(
      wavelength = wavelength,
      stim_nd = stim_nd,
      stim_index = stim_index,
      protocol_id = protocol_id
    ))
  }

  wl <- as.character(wavelength)
  suffix <- if (identical(protocol_id, "C0") && !is.na(stim_index)) {
    if (stim_index >= 1L && stim_index <= 7L) {
      "A"
    } else if (stim_index >= 36L && stim_index <= 42L) {
      "B"
    } else {
      ""
    }
  } else {
    ""
  }

  paste0(wl, suffix, " ", format(as.numeric(stim_nd), nsmall = 1, trim = TRUE))
}

# ==============================================================================
# 1) IMPORT ORIGIN STIMRESP AND BUILD ORIGIN PROTOCOL
# ==============================================================================

stimresp_raw <- readxl::read_excel(
  worksheet_path,
  sheet = stimresp_sheet,
  col_names = TRUE
)

cat("\nOriginal StimResp column names:\n")
print(names(stimresp_raw))

trace_cols <- names(stimresp_raw)[
  stringr::str_detect(
    names(stimresp_raw),
    "^(650\\.00A|650\\.00B|[0-9]+\\.00)\\s+[0-9]+\\.?[0-9]*$"
  )
]

cat("\nNumber of worksheet trace columns:\n")
print(length(trace_cols))

cat("\nWorksheet trace columns in exact order:\n")
print(trace_cols)

origin_protocol_70 <- tibble::tibble(
  stim_index_70 = seq_along(trace_cols),
  stim_col = trace_cols,
  stim_label = canonicalize_stim_label(trace_cols)
) |>
  tidyr::extract(
    stim_label,
    into = c("wavelength", "stim_nd"),
    regex = "^(650A|650B|[0-9]+)\\s+([0-9.]+)$",
    remove = FALSE
  ) |>
  dplyr::mutate(
    wavelength = stringr::str_remove(as.character(wavelength), "[AB]$"),
    stimulus_ND = as.numeric(stim_nd)
  ) |>
  dplyr::select(stim_index_70, stim_col, stim_label, wavelength, stimulus_ND)

cat("\nOrigin-derived 70-label protocol:\n")
print(origin_protocol_70, n = Inf)

origin_protocol_280 <- origin_protocol_70 |>
  tidyr::uncount(weights = 4, .id = "tech_rep") |>
  dplyr::mutate(sweep = dplyr::row_number()) |>
  dplyr::select(sweep, stim_index_70, tech_rep, stim_label, wavelength, stimulus_ND)

cat("\nOrigin-derived 280-sweep protocol:\n")
print(origin_protocol_280, n = Inf)

# ==============================================================================
# 2) STANDARDIZE YOUR UPDATED protocol_table_C0
# ==============================================================================

pt <- protocol_table_C0()

if (is.function(pt)) {
  stop("protocol_table_C0 is a function, not a data frame. Call the constructor first.")
}

cat("\nColumns in protocol_table_C0():\n")
print(names(pt))

get_chr <- function(data, nm) {
  if (nm %in% names(data)) as.character(data[[nm]]) else rep(NA_character_, nrow(data))
}

get_num <- function(data, nm) {
  if (nm %in% names(data)) suppressWarnings(as.numeric(data[[nm]])) else rep(NA_real_, nrow(data))
}

if (!"sweep" %in% names(pt) && !"sweep_in_protocol" %in% names(pt)) {
  pt_sweeps <- expand_protocol_repeats(pt, repeats_per_stim = 4L) |>
    dplyr::mutate(sweep = .data$sweep_in_protocol)
} else if ("sweep_in_protocol" %in% names(pt) && !"sweep" %in% names(pt)) {
  pt_sweeps <- pt |>
    dplyr::mutate(sweep = .data$sweep_in_protocol)
} else {
  pt_sweeps <- pt
}

protocol_new_std <- pt_sweeps |>
  dplyr::mutate(
    sweep = dplyr::coalesce(
      if ("sweep" %in% names(pt_sweeps)) as.integer(pt_sweeps[["sweep"]]) else rep(NA_integer_, nrow(pt_sweeps)),
      dplyr::row_number()
    ),
    stim_index_std = dplyr::coalesce(
      if ("stim_index" %in% names(pt_sweeps)) as.integer(pt_sweeps[["stim_index"]]) else rep(NA_integer_, nrow(pt_sweeps)),
      if ("stim_index_70" %in% names(pt_sweeps)) as.integer(pt_sweeps[["stim_index_70"]]) else rep(NA_integer_, nrow(pt_sweeps))
    ),
    protocol_id_std = dplyr::coalesce(
      get_chr(pt_sweeps, "protocol_id"),
      rep("C0", nrow(pt_sweeps))
    ),
    wavelength_std = dplyr::coalesce(
      get_chr(pt_sweeps, "wavelength"),
      get_chr(pt_sweeps, "wavelength_map"),
      get_chr(pt_sweeps, "wavelength_raw")
    ),
    stim_nd_std = dplyr::coalesce(
      get_num(pt_sweeps, "stim_nd"),
      get_num(pt_sweeps, "stimulus_ND"),
      get_num(pt_sweeps, "stim_nd_map"),
      get_num(pt_sweeps, "stim_nd_raw")
    )
  ) |>
  dplyr::transmute(
    sweep = sweep,
    stim_index_70 = stim_index_std,
    tech_rep = if ("tech_rep" %in% names(pt_sweeps)) as.integer(pt_sweeps[["tech_rep"]]) else NA_integer_,
    wavelength = wavelength_std,
    stimulus_ND = stim_nd_std,
    stim_label = purrr::pmap_chr(
      list(wavelength, stimulus_ND, stim_index_70, protocol_id_std),
      make_stim_label
    )
  ) |>
  dplyr::filter(
    !is.na(sweep),
    !is.na(wavelength),
    !is.na(stimulus_ND)
  ) |>
  dplyr::arrange(sweep)

cat("\nStandardized UPDATED protocol_table_C0 preview:\n")
print(protocol_new_std, n = Inf)

# ==============================================================================
# 3) COMPARE UPDATED PROTOCOL AGAINST ORIGIN PROTOCOL ORDER
# ==============================================================================

protocol_compare <- origin_protocol_280 |>
  dplyr::rename(
    origin_stim_label = stim_label,
    origin_wavelength = wavelength,
    origin_stimulus_ND = stimulus_ND
  ) |>
  dplyr::left_join(
    protocol_new_std |>
      dplyr::rename(
        zeus_stim_label = stim_label,
        zeus_wavelength = wavelength,
        zeus_stimulus_ND = stimulus_ND
      ),
    by = "sweep"
  ) |>
  dplyr::mutate(
    label_match = origin_stim_label == zeus_stim_label,
    wavelength_match = origin_wavelength == zeus_wavelength,
    nd_match = origin_stimulus_ND == zeus_stimulus_ND
  )

cat("\nProtocol comparison summary:\n")
print(
  protocol_compare |>
    dplyr::summarise(
      n_sweeps = dplyr::n(),
      n_label_match = sum(label_match, na.rm = TRUE),
      n_wavelength_match = sum(wavelength_match, na.rm = TRUE),
      n_nd_match = sum(nd_match, na.rm = TRUE)
    )
)

cat("\nProtocol mismatches:\n")
print(
  protocol_compare |>
    dplyr::filter(!(label_match & wavelength_match & nd_match)),
  n = Inf
)

# ==============================================================================
# 4) IMPORT ORIGIN STIMRESP LONG FORMAT
# ==============================================================================

time_col_name <- names(stimresp_raw)[1]

time_diag <- stimresp_raw |>
  dplyr::transmute(
    row_id = dplyr::row_number(),
    time_raw = .data[[time_col_name]],
    time_ms = safe_num(.data[[time_col_name]])
  )

cat("\nTime-column diagnostics:\n")
print(head(time_diag, 25), n = 25)

stimresp_long <- stimresp_raw |>
  dplyr::select(dplyr::all_of(c(time_col_name, trace_cols))) |>
  dplyr::rename(time_raw = dplyr::all_of(time_col_name)) |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    time_ms = safe_num(time_raw)
  ) |>
  dplyr::filter(!is.na(time_ms)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(trace_cols),
    names_to = "stim_col",
    values_to = "worksheet_value"
  ) |>
  dplyr::mutate(
    worksheet_value = safe_num(worksheet_value),
    stim_label = canonicalize_stim_label(stim_col)
  ) |>
  dplyr::left_join(
    origin_protocol_70 |>
      dplyr::select(stim_col, wavelength, stimulus_ND, stim_label),
    by = c("stim_col", "stim_label")
  ) |>
  dplyr::select(stim_label, wavelength, stimulus_ND, time_ms, worksheet_value)

cat("\nParsed StimResp preview:\n")
print(stimresp_long, n = 20)

# ==============================================================================
# 5) IMPORT / STANDARDIZE ZEUS LONG DATA USING UPDATED C0 PROTOCOL
# ==============================================================================

cat("\nZEUS df_long column names:\n")
print(names(df_long))

joined_raw <- df_long |>
  dplyr::filter(channel == channel_keep) |>
  dplyr::left_join(
    protocol_new_std |>
      dplyr::select(sweep, stim_label, wavelength, stimulus_ND),
    by = "sweep"
  )

cat("\nColumns after ZEUS raw join:\n")
print(names(joined_raw))

zeus_trace_long <- joined_raw |>
  dplyr::mutate(
    wavelength_final = dplyr::coalesce(
      safe_chr_col(joined_raw, "wavelength.y"),
      safe_chr_col(joined_raw, "wavelength.x"),
      safe_chr_col(joined_raw, "wavelength")
    ),
    stim_nd_final = dplyr::coalesce(
      safe_num_col(joined_raw, "stimulus_ND"),
      safe_num_col(joined_raw, "stim_nd.y"),
      safe_num_col(joined_raw, "stim_nd.x"),
      safe_num_col(joined_raw, "stim_nd")
    ),
    stim_label_final = dplyr::coalesce(
      safe_chr_col(joined_raw, "stim_label.y"),
      safe_chr_col(joined_raw, "stim_label.x"),
      safe_chr_col(joined_raw, "stim_label")
    )
  ) |>
  dplyr::transmute(
    sweep = sweep,
    stim_label = stim_label_final,
    wavelength = wavelength_final,
    stimulus_ND = stim_nd_final,
    time_ms = safe_num(time_ms),
    value = safe_num(value)
  ) |>
  dplyr::filter(
    !is.na(stim_label),
    !is.na(wavelength),
    !is.na(stimulus_ND),
    !is.na(time_ms),
    !is.na(value)
  )

cat("\nZEUS relabeled trace preview:\n")
print(zeus_trace_long, n = 20)

# ==============================================================================
# 6) BUILD ZEUS MEAN TRACE FROM UPDATED PROTOCOL
# ==============================================================================

zeus_mean_trace <- zeus_trace_long |>
  dplyr::group_by(stim_label, wavelength, stimulus_ND, time_ms) |>
  dplyr::summarise(
    zeus_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nWorksheet-aligned ZEUS mean trace preview:\n")
print(zeus_mean_trace, n = 20)

# ==============================================================================
# 7) JOIN ZEUS MEAN TRACE TO ORIGIN STIMRESP TRACE
# ==============================================================================

wave_compare <- zeus_mean_trace |>
  dplyr::inner_join(
    stimresp_long,
    by = c("stim_label", "wavelength", "stimulus_ND", "time_ms")
  ) |>
  dplyr::mutate(
    diff = zeus_value - worksheet_value,
    abs_diff = abs(diff),
    sq_diff = diff^2
  )

cat("\nRows in waveform comparison:\n")
print(nrow(wave_compare))

if (nrow(wave_compare) == 0) {
  stop("No matched rows in waveform comparison. Check protocol labeling and time alignment.", call. = FALSE)
}

# ==============================================================================
# 8) TRACE AGREEMENT SUMMARY
# ==============================================================================

overall_wave_summary <- wave_compare |>
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor_value = suppressWarnings(cor(zeus_value, worksheet_value, use = "complete.obs"))
  )

cat("\nOverall waveform agreement:\n")
print(overall_wave_summary)

by_stim_wave_summary <- wave_compare |>
  dplyr::group_by(stim_label, wavelength, stimulus_ND) |>
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor_value = suppressWarnings(cor(zeus_value, worksheet_value, use = "complete.obs")),
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_abs_diff)

cat("\nBy-stimulus waveform agreement:\n")
print(by_stim_wave_summary, n = Inf)

worst_labels <- by_stim_wave_summary |>
  dplyr::arrange(dplyr::desc(mean_abs_diff))

cat("\nWorst labels:\n")
print(worst_labels, n = 20)

worst_points <- wave_compare |>
  dplyr::arrange(dplyr::desc(abs_diff)) |>
  dplyr::slice_head(n = 25)

cat("\nWorst waveform points:\n")
print(worst_points, n = 25)

# ==============================================================================
# 9) OPTIONAL: BASELINE-ZEROED TRACE COMPARISON
# Useful if raw-vs-Origin mismatch is mostly baseline offset
# ==============================================================================

baseline_zero_by_label <- function(df, value_col) {
  df |>
    dplyr::group_by(stim_label, wavelength, stimulus_ND) |>
    dplyr::mutate(
      baseline_est = mean(
        .data[[value_col]][.data$time_ms >= 300 & .data$time_ms <= 400],
        na.rm = TRUE
      ),
      value_zeroed = .data[[value_col]] - baseline_est
    ) |>
    dplyr::ungroup()
}

zeus_zeroed <- baseline_zero_by_label(
  zeus_mean_trace |>
    dplyr::rename(value = zeus_value),
  "value"
) |>
  dplyr::rename(zeus_zeroed = value_zeroed) |>
  dplyr::select(stim_label, wavelength, stimulus_ND, time_ms, zeus_zeroed)

origin_zeroed <- baseline_zero_by_label(
  stimresp_long |>
    dplyr::rename(value = worksheet_value),
  "value"
) |>
  dplyr::rename(origin_zeroed = value_zeroed) |>
  dplyr::select(stim_label, wavelength, stimulus_ND, time_ms, origin_zeroed)

wave_compare_zeroed <- zeus_zeroed |>
  dplyr::inner_join(
    origin_zeroed,
    by = c("stim_label", "wavelength", "stimulus_ND", "time_ms")
  ) |>
  dplyr::mutate(
    diff = zeus_zeroed - origin_zeroed,
    abs_diff = abs(diff),
    sq_diff = diff^2
  )

baseline_zeroed_summary <- wave_compare_zeroed |>
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor_value = suppressWarnings(cor(zeus_zeroed, origin_zeroed, use = "complete.obs"))
  )

cat("\nBaseline-zeroed waveform agreement:\n")
print(baseline_zeroed_summary)

baseline_zeroed_by_label <- wave_compare_zeroed |>
  dplyr::group_by(stim_label) |>
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor_value = suppressWarnings(cor(zeus_zeroed, origin_zeroed, use = "complete.obs")),
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_abs_diff)

cat("\nBaseline-zeroed waveform agreement by stim label:\n")
print(baseline_zeroed_by_label, n = Inf)

# ==============================================================================
# 10) OBJECTS TO KEEP
# ==============================================================================

comparison_results <- list(
  origin_protocol_70 = origin_protocol_70,
  origin_protocol_280 = origin_protocol_280,
  protocol_new_std = protocol_new_std,
  protocol_compare = protocol_compare,
  stimresp_long = stimresp_long,
  zeus_trace_long = zeus_trace_long,
  zeus_mean_trace = zeus_mean_trace,
  wave_compare = wave_compare,
  overall_wave_summary = overall_wave_summary,
  by_stim_wave_summary = by_stim_wave_summary,
  wave_compare_zeroed = wave_compare_zeroed,
  baseline_zeroed_summary = baseline_zeroed_summary,
  baseline_zeroed_by_label = baseline_zeroed_by_label,
  worst_points = worst_points
)

cat("\nFinished. Results stored in `comparison_results`.\n")
