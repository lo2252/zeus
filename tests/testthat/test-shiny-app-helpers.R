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
