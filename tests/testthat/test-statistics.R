test_that("zeus_avg_peak_by_nd summarizes peaks across ND", {
  x <- tibble::tibble(
    protocol_id = c("C0", "C0", "C0", "C1"),
    wavelength = c("650", "570", "650", "White"),
    stim_nd = c(4.0, 4.0, 3.0, 6.0),
    amp_mv = c(10, 14, 6, 8),
    awave_mv = c(-4, -6, -2, -3),
    dwave_mv = c(3, 5, 2, 4),
    trough_time_poststim_ms = c(20, 22, 18, 21),
    peak_time_poststim_ms = c(50, 54, 48, 52),
    dpeak_time_poststim_ms = c(120, 122, 118, 121)
  )

  out <- zeus_avg_peak_by_nd(x)

  expect_true(all(c("protocol_id", "stim_nd", "peak_type", "mean_peak_mv", "sd_peak_mv") %in% names(out)))

  c0_bwave_nd4 <- out |>
    dplyr::filter(.data$protocol_id == "C0", .data$stim_nd == 4, .data$peak_type == "B-wave")

  expect_equal(c0_bwave_nd4$n, 2)
  expect_equal(c0_bwave_nd4$mean_peak_mv, 12)
  expect_equal(round(c0_bwave_nd4$sd_peak_mv, 6), round(stats::sd(c(10, 14)), 6))
})

test_that("zeus_avg_peak_by_wavelength excludes white-light summaries", {
  x <- tibble::tibble(
    protocol_id = c("C0", "C0", "C0", "C1"),
    wavelength = c("650", "650", "570", "White"),
    stim_label = c("650A 4.0", "650B 3.0", "570 4.0", "White 6.0"),
    stim_nd = c(4.0, 3.0, 4.0, 6.0),
    amp_mv = c(10, 6, 14, 8),
    awave_mv = c(-4, -2, -6, -3),
    dwave_mv = c(3, 2, 5, 4),
    trough_time_poststim_ms = c(20, 18, 22, 21),
    peak_time_poststim_ms = c(50, 48, 54, 52),
    dpeak_time_poststim_ms = c(120, 118, 122, 121)
  )

  out <- zeus_avg_peak_by_wavelength(x)

  expect_false(any(out$wavelength == "White"))
  expect_true(all(c("650A", "650B", "570") %in% out$wavelength))

  c0_650a_bwave <- out |>
    dplyr::filter(.data$protocol_id == "C0", .data$wavelength == "650A", .data$peak_type == "B-wave")

  expect_equal(c0_650a_bwave$n, 1)
  expect_equal(c0_650a_bwave$mean_peak_mv, 10)
})

test_that("zeus_summarize_peak_statistics returns report-ready components", {
  x <- tibble::tibble(
    protocol_id = c("C0", "C0", "C1"),
    wavelength = c("650", "570", "White"),
    stim_nd = c(4.0, 4.0, 6.0),
    amp_mv = c(10, 14, 8),
    awave_mv = c(-4, -6, -3),
    dwave_mv = c(3, 5, 4),
    trough_time_poststim_ms = c(20, 22, 21),
    peak_time_poststim_ms = c(50, 54, 52),
    dpeak_time_poststim_ms = c(120, 122, 121)
  )

  out <- zeus_summarize_peak_statistics(x)

  expect_named(out, c("key_statistics", "by_nd", "by_wavelength", "combined_export"))
  expect_true(nrow(out$key_statistics) > 0)
  expect_true(nrow(out$by_nd) > 0)
  expect_true(nrow(out$by_wavelength) > 0)
  expect_true(nrow(out$combined_export) > 0)
  expect_true(all(c("key_statistics", "by_nd", "by_wavelength") %in% out$combined_export$summary_type))
})

test_that("zeus_summarize_peak_statistics handles all-missing peak groups without warnings", {
  x <- tibble::tibble(
    protocol_id = c("C1", "C1"),
    wavelength = c("White", "White"),
    stim_nd = c(6.0, 5.5),
    amp_mv = c(8, 7),
    awave_mv = c(-3, -2),
    dwave_mv = c(NA_real_, NA_real_),
    trough_time_poststim_ms = c(21, 22),
    peak_time_poststim_ms = c(52, 53),
    dpeak_time_poststim_ms = c(NA_real_, NA_real_)
  )

  expect_no_warning(out <- zeus_summarize_peak_statistics(x))
  expect_true(all(is.na(out$key_statistics$min_peak_mv[out$key_statistics$peak_type == "D-wave"])))
})
