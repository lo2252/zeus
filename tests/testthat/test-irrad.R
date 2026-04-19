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

test_that("waveform marker helper uses dedicated A-wave timing", {
  time_ms <- seq(300, 1000, by = 1)
  overall_df <- tibble::tibble(
    time_ms = time_ms,
    signal = 0.03 * sin((time_ms - 400) / 20) -
      0.55 * exp(-((time_ms - 455) / 12)^2) +
      1.2 * exp(-((time_ms - 560) / 24)^2) +
      0.4 * exp(-((time_ms - 860) / 25)^2)
  )

  markers <- .zeus_waveform_marker_lines(
    overall_df = overall_df,
    a_window = c(430, 500),
    b_window = c(400, 700),
    d_window = c(700, 1000)
  )

  expect_equal(as.character(markers$marker_type), c("A-wave trough", "B-wave peak", "D-wave peak"))
  expect_equal(markers$time_ms[[1]], 455)
  expect_true(markers$time_ms[[2]] > markers$time_ms[[1]])
})

test_that("waveform marker helper keeps A-wave before the B-wave peak", {
  time_ms <- seq(300, 1000, by = 1)
  overall_df <- tibble::tibble(
    time_ms = time_ms,
    signal = 0.02 * sin((time_ms - 400) / 18) -
      0.20 * exp(-((time_ms - 450) / 10)^2) +
      1.10 * exp(-((time_ms - 520) / 18)^2) -
      0.80 * exp(-((time_ms - 640) / 16)^2) +
      0.35 * exp(-((time_ms - 860) / 25)^2)
  )

  markers <- .zeus_waveform_marker_lines(
    overall_df = overall_df,
    a_window = c(400, 700),
    b_window = c(400, 700),
    d_window = c(700, 1000)
  )

  expect_equal(markers$time_ms[[1]], 450)
  expect_true(markers$time_ms[[1]] < markers$time_ms[[2]])
})
