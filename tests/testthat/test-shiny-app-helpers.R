test_that("app wavelength helper returns non-white C0 blocks", {
  x <- list(
    traces_70 = data.frame(
      stim_label = c("650A 4.0", "570 5.0", "White 6.0", "650A 3.5"),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(
    .zeus_app_available_wavelengths(x),
    c("650A", "570")
  )
})

test_that("app export path helper creates expected file names", {
  out <- .zeus_app_export_paths("results", "fish01")

  expect_equal(out$mean_plot, file.path("results", "fish01_mean_waveform.png"))
  expect_equal(out$spectral_plot, file.path("results", "fish01_spectral_waveform.png"))
  expect_equal(out$intensity_plot, file.path("results", "fish01_intensity_response.png"))
  expect_equal(out$csv_base, file.path("results", "fish01"))
})

test_that("peak settings helper includes awave window", {
  input <- list(
    peak_baseline_start = 300,
    peak_baseline_end = 400,
    peak_response_start = 400,
    peak_response_end = 700,
    awave_window_start = 400,
    awave_window_end = 700,
    stimulus_onset_ms = 400,
    same_sign = TRUE,
    time_reference = "absolute",
    dwave_window_start = 700,
    dwave_window_end = 1000
  )

  out <- .zeus_app_peak_settings(input)

  expect_equal(out$awave_window_ms, c(400, 700))
})

test_that("full fish protocol defaults and labels are assigned predictably", {
  expect_equal(
    .zeus_app_default_protocols(3, analysis_mode = "full_fish"),
    c("C1", "C0", "C0")
  )

  expect_equal(
    .zeus_app_make_item_labels(c("C1", "C0", "C0")),
    c("C1", "C0 A", "C0 B")
  )
})

test_that("full fish protocol validation enforces one C1 and two C0 files", {
  expect_null(.zeus_app_validate_protocols(c("C1", "C0", "C0"), analysis_mode = "full_fish"))
  expect_match(
    .zeus_app_validate_protocols(c("C1", "C1", "C0"), analysis_mode = "full_fish"),
    "1 C1 file and 2 C0 files"
  )
  expect_match(
    .zeus_app_validate_protocols(c("C1", "C0"), analysis_mode = "full_fish"),
    "exactly 3 ABF files"
  )
})

test_that("zeus_read_abf preserves raw traces in traces_280", {
  path <- zeus_extdata("26225004.abf")
  skip_if_not(nzchar(path) && file.exists(path), "26225004.abf not found in inst/extdata")

  x <- zeus_read_abf(path, protocol = "C1", smooth_n = 3L)

  expect_true("value_raw" %in% names(x$traces_280))
  expect_true(any(is.finite(x$traces_280$value_raw)))
  expect_true(any(abs(x$traces_280$value - x$traces_280$value_raw) > 0, na.rm = TRUE))
})
