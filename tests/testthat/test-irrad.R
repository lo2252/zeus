test_that("measure_trace_window returns one row", {
  trace_df <- tibble::tibble(
    time = seq(0, 1, length.out = 100),
    value = c(rep(0, 40), seq(0, -5, length.out = 20), seq(-5, 3, length.out = 40))
  )
  
  out <- measure_trace_window(
    trace_df,
    baseline_window_ms = c(300, 400),
    response_window_ms = c(400, 700),
    stimulus_onset_ms = 400
  )
  
  expect_equal(nrow(out), 1L)
  expect_true(all(c("amp_mv", "trough_time_ms", "peak_time_ms") %in% names(out)))
})

test_that("nd_to_log_irradiance joins calibration values", {
  meta <- tibble::tibble(
    wavelength = c("White", "650"),
    stim_nd = c(6.0, 4.0)
  )
  
  calib <- tibble::tibble(
    wavelength = c("White", "650"),
    stim_nd = c(6.0, 4.0),
    log_hv = c(-5.0, -4.2)
  )
  
  out <- nd_to_log_irradiance(meta, calib)
  expect_equal(out$log_hv, c(-5.0, -4.2))
})