# CSV export tests --------------------------------------------------------

test_that("zeus_export_csv_bundle writes one csv per supported component", {
  obj <- list(
    raw = data.frame(
      sweep = c(1L, 1L),
      time = c(0, 0.001),
      channel = c("ERG DAM80", "ERG DAM80"),
      units = c("uV", "uV"),
      value = c(1.2, 1.5)
    ),
    traces_280 = data.frame(
      sweep = c(1L, 2L),
      stim_index = c(1L, 1L),
      value = c(2.1, 2.2)
    ),
    traces_70 = data.frame(
      stim_index = 1L,
      stim_label = "White 3.0",
      time_ms = 100,
      value = 2.15
    ),
    photocell = data.frame(
      sweep = 1L,
      time = 0.001,
      value = 0.5
    ),
    stimresp_qc = data.frame(
      sweep = c(1L, 2L),
      keep = c(TRUE, FALSE)
    ),
    stimresp_settings = list(
      stimresp_zero_baseline = TRUE,
      stimresp_baseline_window_ms = c(300, 400),
      stimresp_runmean_k = 16L
    )
  )

  csv_path <- tempfile(fileext = ".csv")
  base_path <- sub("\\.csv$", "", csv_path, ignore.case = TRUE)

  written <- zeus_export_csv_bundle(obj, csv_path)

  expected <- c(
    raw = paste0(base_path, "_raw.csv"),
    traces_280 = paste0(base_path, "_traces_280.csv"),
    traces_70 = paste0(base_path, "_traces_70.csv"),
    photocell = paste0(base_path, "_photocell.csv"),
    stimresp_qc = paste0(base_path, "_stimresp_qc.csv"),
    stimresp_settings = paste0(base_path, "_stimresp_settings.csv")
  )

  expect_equal(unname(written), unname(expected))
  expect_true(all(file.exists(expected)))

  settings_export <- utils::read.csv(
    expected[["stimresp_settings"]],
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(settings_export), 1L)
  expect_true(all(c(
    "stimresp_zero_baseline",
    "stimresp_baseline_window_ms_1",
    "stimresp_baseline_window_ms_2",
    "stimresp_runmean_k"
  ) %in% names(settings_export)))
})