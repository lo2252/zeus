test_that("make_stimresp_label adds A/B suffixes for C0", {
  expect_equal(make_stimresp_label("650", 4.0, stim_index = 1, protocol_id = "C0"), "650A 4.0")
  expect_equal(make_stimresp_label("650", 4.5, stim_index = 36, protocol_id = "C0"), "650B 4.5")
  expect_equal(make_stimresp_label("White", 6.0, stim_index = 1, protocol_id = "C1"), "White 6.0")
})

test_that("zero_baseline_traces preserves row count", {
  df <- tibble::tibble(
    sweep = rep(1:2, each = 5),
    time = rep(c(0.30, 0.35, 0.40, 0.45, 0.50), 2),
    channel = "ERG DAM80",
    value = c(1, 1, 1, 2, 3, 2, 2, 2, 4, 5)
  )
  
  out <- zero_baseline_traces(df, baseline_window_ms = c(300, 400))
  expect_equal(nrow(out), nrow(df))
})