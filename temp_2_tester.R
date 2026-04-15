# ==============================================================================
# ZEUS C1 VALIDATION AGAINST ORIGIN
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(tibble)

if (!exists("zeus_read_abf", mode = "function") ||
    !exists("protocol_table_C1", mode = "function")) {
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

# ------------------------------------------------------------------------------
# INPUTS
# ------------------------------------------------------------------------------

worksheet_path <- "../test_data/Origin_Files_with_dwave/26225004_origin_export_with_d_wave.xlsx"
stimresp_sheet <- "StimResp"
channel_keep   <- "ERG DAM80"
zeus_obj <- zeus_read_abf(
  "../test_data/Origin_Files_with_dwave/26225004.abf",
  protocol = "C1"
)
df_long <- zeus_obj$traces_280
# ------------------------------------------------------------------------------
# CHECK OBJECTS
# ------------------------------------------------------------------------------

if (!exists("protocol_table_C1", mode = "function")) {
  stop("protocol_table_C1 not found.", call. = FALSE)
}

if (!exists("df_long")) {
  stop("df_long not found.", call. = FALSE)
}

stopifnot(is.data.frame(df_long))

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_chr_col <- function(data, nm) {
  if (nm %in% names(data)) as.character(data[[nm]]) else rep(NA_character_, nrow(data))
}

safe_num_col <- function(data, nm) {
  if (nm %in% names(data)) suppressWarnings(as.numeric(data[[nm]])) else rep(NA_real_, nrow(data))
}

canonicalize_stim_label <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace("\\.\\.\\.[0-9]+$", "") |>
    stringr::str_replace("^([0-9]+)\\.00", "\\1") |>
    stringr::str_squish()
}

make_stim_label <- function(wavelength, stim_nd, stim_index = NA_integer_, protocol_id = "C1") {
  if (exists("make_stimresp_label", mode = "function")) {
    return(make_stimresp_label(
      wavelength = wavelength,
      stim_nd = stim_nd,
      stim_index = stim_index,
      protocol_id = protocol_id
    ))
  }

  paste0(as.character(wavelength), " ", format(as.numeric(stim_nd), nsmall = 1, trim = TRUE))
}

# ------------------------------------------------------------------------------
# IMPORT ORIGIN
# ------------------------------------------------------------------------------

stimresp_raw <- readxl::read_excel(
  worksheet_path,
  sheet = stimresp_sheet
)

time_col_name <- names(stimresp_raw)[1]
pre_mean_cols <- names(stimresp_raw)[seq_len(match("mean", names(stimresp_raw)) - 1L)]
trace_cols <- pre_mean_cols[
  stringr::str_detect(pre_mean_cols, "^White\\s+[0-9]+\\.?[0-9]*(\\.\\.\\.[0-9]+)?$")
]

cat("\nNumber of worksheet C1 trace columns:\n")
print(length(trace_cols))

cat("\nWorksheet C1 trace columns in exact order:\n")
print(trace_cols)

# ------------------------------------------------------------------------------
# BUILD ORIGIN PROTOCOL (C1 = SAME WAVELENGTH)
# ------------------------------------------------------------------------------

origin_protocol_70 <- tibble(
  stim_index_70 = seq_along(trace_cols),
  stim_col = trace_cols,
  stim_label = canonicalize_stim_label(trace_cols)
) |>
  tidyr::extract(
    stim_label,
    into = c("wavelength", "stim_nd"),
    regex = "^(White)\\s+([0-9.]+)$",
    remove = FALSE
  ) |>
  mutate(
    wavelength = as.character(wavelength),
    stimulus_ND = as.numeric(stim_nd)
  ) |>
  select(stim_index_70, stim_col, stim_label, wavelength, stimulus_ND)

origin_protocol_280 <- origin_protocol_70 |>
  tidyr::uncount(weights = 4, .id = "tech_rep") |>
  mutate(sweep = row_number()) |>
  select(sweep, stim_index_70, tech_rep, stim_label, wavelength, stimulus_ND)

# ------------------------------------------------------------------------------
# STANDARDIZE C1 PROTOCOL
# ------------------------------------------------------------------------------

pt <- protocol_table_C1()

pt_sweeps <- expand_protocol_repeats(pt, repeats_per_stim = 4L) |>
  mutate(sweep = .data$sweep_in_protocol)

protocol_new_std <- pt_sweeps |>
  mutate(
    stim_index_70 = as.integer(stim_index),
    wavelength = as.character(wavelength),
    stimulus_ND = as.numeric(stim_nd),
    stim_label = purrr::pmap_chr(
      list(wavelength, stimulus_ND, stim_index_70, protocol_id),
      make_stim_label
    )
  ) |>
  select(sweep, stim_index_70, tech_rep, stim_label, wavelength, stimulus_ND)

# ------------------------------------------------------------------------------
# PROTOCOL CHECK
# ------------------------------------------------------------------------------

protocol_compare <- origin_protocol_280 |>
  rename(
    origin_stim_index_70 = stim_index_70,
    origin_tech_rep = tech_rep,
    origin_stim_label = stim_label,
    origin_wavelength = wavelength,
    origin_stimulus_ND = stimulus_ND
  ) |>
  left_join(
    protocol_new_std |>
      rename(
        zeus_stim_index_70 = stim_index_70,
        zeus_tech_rep = tech_rep,
        zeus_stim_label = stim_label,
        zeus_wavelength = wavelength,
        zeus_stimulus_ND = stimulus_ND
      ),
    by = "sweep"
  ) |>
  mutate(
    index_match = origin_stim_index_70 == zeus_stim_index_70,
    label_match = origin_stim_label == zeus_stim_label,
    wavelength_match = origin_wavelength == zeus_wavelength,
    nd_match = origin_stimulus_ND == zeus_stimulus_ND
  )

cat("\nC1 Protocol match summary:\n")
print(
  protocol_compare |>
    summarise(
      n_match = sum(index_match & label_match & wavelength_match & nd_match, na.rm = TRUE),
      n_total = n()
    )
)

cat("\nC1 protocol mismatches:\n")
print(
  protocol_compare |>
    filter(!(index_match & label_match & wavelength_match & nd_match)),
  n = Inf
)

# ------------------------------------------------------------------------------
# BUILD ORIGIN LONG
# ------------------------------------------------------------------------------

stimresp_long <- stimresp_raw |>
  select(all_of(c(time_col_name, trace_cols))) |>
  rename(time_raw = all_of(time_col_name)) |>
  mutate(time_ms = safe_num(time_raw)) |>
  filter(!is.na(time_ms)) |>
  pivot_longer(
    cols = all_of(trace_cols),
    names_to = "stim_col",
    values_to = "origin_value"
  ) |>
  mutate(
    origin_value = safe_num(origin_value),
    stim_label = canonicalize_stim_label(stim_col)
  ) |>
  left_join(
    origin_protocol_70 |>
      select(stim_col, stim_index_70, wavelength, stimulus_ND, stim_label),
    by = c("stim_col", "stim_label")
  ) |>
  select(stim_index_70, stim_label, wavelength, stimulus_ND, time_ms, origin_value)

# ------------------------------------------------------------------------------
# BUILD ZEUS TRACE
# ------------------------------------------------------------------------------

joined <- df_long |>
  filter(channel == channel_keep) |>
  left_join(
    protocol_new_std |>
      select(sweep, stim_index_70, stim_label, wavelength, stimulus_ND),
    by = "sweep"
  )

zeus_trace <- joined |>
  mutate(
    stim_index_70_final = dplyr::coalesce(
      safe_num_col(joined, "stim_index_70"),
      safe_num_col(joined, "stim_index")
    ),
    stim_label_final = dplyr::coalesce(
      safe_chr_col(joined, "stim_label.y"),
      safe_chr_col(joined, "stim_label.x"),
      safe_chr_col(joined, "stim_label")
    ),
    wavelength_final = dplyr::coalesce(
      safe_chr_col(joined, "wavelength.y"),
      safe_chr_col(joined, "wavelength.x"),
      safe_chr_col(joined, "wavelength")
    ),
    stim_nd_final = dplyr::coalesce(
      safe_num_col(joined, "stimulus_ND"),
      safe_num_col(joined, "stim_nd.y"),
      safe_num_col(joined, "stim_nd.x"),
      safe_num_col(joined, "stim_nd")
    )
  ) |>
  transmute(
    sweep = sweep,
    stim_index_70 = as.integer(stim_index_70_final),
    stim_label = stim_label_final,
    wavelength = wavelength_final,
    stimulus_ND = stim_nd_final,
    time_ms = safe_num(time_ms),
    value = safe_num(value)
  ) |>
  filter(
    !is.na(stim_index_70),
    !is.na(stim_label),
    !is.na(wavelength),
    !is.na(stimulus_ND),
    !is.na(time_ms),
    !is.na(value)
  )

zeus_mean <- zeus_trace |>
  group_by(stim_index_70, stim_label, wavelength, stimulus_ND, time_ms) |>
  summarise(zeus_value = mean(value, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# COMPARE
# ------------------------------------------------------------------------------

wave_compare <- zeus_mean |>
  inner_join(
    stimresp_long,
    by = c("stim_index_70", "stim_label", "wavelength", "stimulus_ND", "time_ms")
  ) |>
  mutate(
    diff = zeus_value - origin_value,
    abs_diff = abs(diff),
    sq_diff = diff^2
  )

summary_results <- wave_compare |>
  summarise(
    n_points = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor = cor(zeus_value, origin_value, use = "complete.obs")
  )

cat("\nC1 waveform agreement:\n")
print(summary_results)

# ------------------------------------------------------------------------------
# WORST STIMS
# ------------------------------------------------------------------------------

by_stim <- wave_compare |>
  group_by(stim_index_70, stim_label, wavelength, stimulus_ND) |>
  summarise(
    n_points = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor = cor(zeus_value, origin_value, use = "complete.obs"),
    .groups = "drop"
  ) |>
  arrange(desc(mean_abs_diff))

cat("\nWorst C1 stim labels:\n")
print(by_stim, n = 10)

# ------------------------------------------------------------------------------
# OPTIONAL: BASELINE-ZEROED COMPARISON
# ------------------------------------------------------------------------------

baseline_zero_by_stim <- function(df, value_col) {
  df |>
    group_by(stim_index_70, stim_label, wavelength, stimulus_ND) |>
    mutate(
      baseline_est = mean(
        .data[[value_col]][.data$time_ms >= 300 & .data$time_ms <= 400],
        na.rm = TRUE
      ),
      value_zeroed = .data[[value_col]] - baseline_est
    ) |>
    ungroup()
}

zeus_zeroed <- baseline_zero_by_stim(
  zeus_mean |> rename(value = zeus_value),
  "value"
) |>
  rename(zeus_zeroed = value_zeroed) |>
  select(stim_index_70, stim_label, wavelength, stimulus_ND, time_ms, zeus_zeroed)

origin_zeroed <- baseline_zero_by_stim(
  stimresp_long |> rename(value = origin_value),
  "value"
) |>
  rename(origin_zeroed = value_zeroed) |>
  select(stim_index_70, stim_label, wavelength, stimulus_ND, time_ms, origin_zeroed)

wave_compare_zeroed <- zeus_zeroed |>
  inner_join(
    origin_zeroed,
    by = c("stim_index_70", "stim_label", "wavelength", "stimulus_ND", "time_ms")
  ) |>
  mutate(
    diff = zeus_zeroed - origin_zeroed,
    abs_diff = abs(diff),
    sq_diff = diff^2
  )

baseline_zeroed_summary <- wave_compare_zeroed |>
  summarise(
    n_points = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(sq_diff, na.rm = TRUE)),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    cor = cor(zeus_zeroed, origin_zeroed, use = "complete.obs")
  )

cat("\nC1 baseline-zeroed waveform agreement:\n")
print(baseline_zeroed_summary)

comparison_results <- list(
  origin_protocol_70 = origin_protocol_70,
  origin_protocol_280 = origin_protocol_280,
  protocol_new_std = protocol_new_std,
  protocol_compare = protocol_compare,
  stimresp_long = stimresp_long,
  zeus_trace = zeus_trace,
  zeus_mean = zeus_mean,
  wave_compare = wave_compare,
  summary_results = summary_results,
  by_stim = by_stim,
  wave_compare_zeroed = wave_compare_zeroed,
  baseline_zeroed_summary = baseline_zeroed_summary
)

cat("\nFinished. Results stored in `comparison_results`.\n")
