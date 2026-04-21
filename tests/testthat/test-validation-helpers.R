# Validation helper tests -------------------------------------------------

test_that("zeus_summarize_protocol_validation summarizes logical match columns", {
  compare_df <- data.frame(
    label_match = c(TRUE, TRUE, FALSE, TRUE),
    wavelength_match = c(TRUE, TRUE, TRUE, TRUE),
    nd_match = c(TRUE, FALSE, TRUE, TRUE)
  )

  out <- zeus_summarize_protocol_validation(compare_df)

  expect_equal(nrow(out), 3L)
  expect_equal(out$n_match[out$comparison == "label_match"], 3)
  expect_equal(out$n_total[out$comparison == "label_match"], 4)
  expect_equal(out$pct_match[out$comparison == "wavelength_match"], 100)
})

test_that("zeus_validate_protocol_agreement requires complete matches", {
  compare_df <- data.frame(match = c(TRUE, TRUE, TRUE))

  out <- zeus_validate_protocol_agreement(compare_df, match_cols = "match")

  expect_true(out$passed)
  expect_error(
    zeus_validate_protocol_agreement(
      data.frame(match = c(TRUE, FALSE)),
      match_cols = "match"
    ),
    "Protocol validation failed"
  )
})

test_that("zeus_compare_waveforms returns overall and grouped metrics", {
  zeus_df <- data.frame(
    stim_label = c("A", "A", "B", "B"),
    time_ms = c(0, 1, 0, 1),
    zeus_value = c(1, 2, 2, 4)
  )

  ref_df <- data.frame(
    stim_label = c("A", "A", "B", "B"),
    time_ms = c(0, 1, 0, 1),
    reference_value = c(1, 1, 3, 5)
  )

  out <- zeus_compare_waveforms(
    zeus_df = zeus_df,
    reference_df = ref_df
  )

  expect_true(all(c(
    "point_compare",
    "overall_summary",
    "group_summary"
  ) %in% names(out)))

  expect_equal(nrow(out$point_compare), 4L)
  expect_equal(out$overall_summary$n_points, 4L)
  expect_equal(round(out$overall_summary$mean_abs_diff, 6), 0.75)
  expect_equal(nrow(out$group_summary), 2L)
  expect_true(all(c("A", "B") %in% out$group_summary$stim_label))
  expect_equal(out$coverage_summary$n_matched_rows, 4L)
  expect_equal(out$coverage_summary$pct_zeus_matched, 100)
  expect_true("n_complete" %in% names(out$overall_summary))
})

test_that("zeus_compare_waveforms rejects duplicate keys", {
  zeus_df <- data.frame(
    stim_label = c("A", "A"),
    time_ms = c(0, 0),
    zeus_value = c(1, 2)
  )

  ref_df <- data.frame(
    stim_label = "A",
    time_ms = 0,
    reference_value = 1
  )

  expect_error(
    zeus_compare_waveforms(zeus_df, ref_df),
    "duplicate join keys"
  )
})

test_that("zeus_compare_waveforms requires numeric response values", {
  zeus_df <- data.frame(
    stim_label = "A",
    time_ms = 0,
    zeus_value = "1"
  )

  ref_df <- data.frame(
    stim_label = "A",
    time_ms = 0,
    reference_value = 1
  )

  expect_error(
    zeus_compare_waveforms(zeus_df, ref_df),
    "must be numeric"
  )
})

test_that("zeus_validate_response_agreement reports threshold checks", {
  summary_df <- data.frame(
    n_points = 4L,
    n_complete = 4L,
    mean_abs_diff = 0.75,
    rmse = sqrt(0.75),
    max_abs_diff = 1,
    cor_value = 0.9
  )

  out <- zeus_validate_response_agreement(
    summary_df,
    min_n_complete = 4L,
    min_correlation = 0.8,
    max_mean_abs_diff = 1,
    max_rmse = 1,
    max_abs_diff = 2,
    label = "test response"
  )

  expect_equal(nrow(out), 5L)
  expect_true(all(out$passed))

  expect_error(
    zeus_validate_response_agreement(
      summary_df,
      min_correlation = 0.95,
      label = "test response"
    ),
    "Response validation failed"
  )

  summary_df$cor_value <- NA_real_

  expect_error(
    zeus_validate_response_agreement(
      summary_df,
      min_correlation = 0.8,
      label = "test response"
    ),
    "Response validation failed"
  )
})
