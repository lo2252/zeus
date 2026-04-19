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
})
