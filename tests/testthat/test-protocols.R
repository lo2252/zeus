test_that("protocol_table_C0 returns 70 rows", {
  x <- protocol_table_C0()
  
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), 70L)
  expect_equal(length(unique(x$block_index)), 10L)
  expect_equal(table(x$block_index)[1], 7L)
})

test_that("protocol_table_C1 returns 70 rows", {
  x <- protocol_table_C1()
  
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), 70L)
  expect_true(all(x$wavelength == "White"))
})

test_that("expand_protocol_repeats returns 280 rows with 4 repeats", {
  x <- expand_protocol_repeats(protocol_table_C0(), repeats_per_stim = 4L)
  
  expect_equal(nrow(x), 280L)
  expect_equal(max(x$tech_rep), 4L)
})