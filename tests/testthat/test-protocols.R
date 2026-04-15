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

# Origin-validated C0 protocol label tests --------------------------------
# These labels were exported directly from Origin StimResp and are used to
# confirm that the Zeus protocol transformation produces the same 70
# stimulus conditions in the same order.

test_that("C0 protocol produces stim_labels matching Origin export", {
  # Expected 70 stim_labels in Origin column order (validated against
  # temp_file/comparison_results_export/origin_protocol_70.csv)
  origin_labels <- c(
    "650A 4.0", "650A 3.5", "650A 3.0", "650A 2.5", "650A 2.0", "650A 1.5", "650A 1.0",
    "570 5.0",  "570 4.5",  "570 4.0",  "570 3.5",  "570 3.0",  "570 2.5",  "570 2.0",
    "490 5.0",  "490 4.5",  "490 4.0",  "490 3.5",  "490 3.0",  "490 2.5",  "490 2.0",
    "410 5.0",  "410 4.5",  "410 4.0",  "410 3.5",  "410 3.0",  "410 2.5",  "410 2.0",
    "330 4.0",  "330 3.5",  "330 3.0",  "330 2.5",  "330 2.0",  "330 1.5",  "330 1.0",
    "650B 4.5", "650B 4.0", "650B 3.5", "650B 3.0", "650B 2.5", "650B 2.0", "650B 1.5",
    "610 5.0",  "610 4.5",  "610 4.0",  "610 3.5",  "610 3.0",  "610 2.5",  "610 2.0",
    "530 5.0",  "530 4.5",  "530 4.0",  "530 3.5",  "530 3.0",  "530 2.5",  "530 2.0",
    "450 5.0",  "450 4.5",  "450 4.0",  "450 3.5",  "450 3.0",  "450 2.5",  "450 2.0",
    "370 5.0",  "370 4.5",  "370 4.0",  "370 3.5",  "370 3.0",  "370 2.5",  "370 2.0"
  )

  tbl <- protocol_table_C0()

  zeus_labels <- purrr::pmap_chr(
    list(tbl$wavelength, tbl$stim_nd, tbl$stim_index, tbl$protocol_id),
    make_stimresp_label
  )

  expect_equal(zeus_labels, origin_labels)
})

test_that("C0 650 nm blocks have correct A/B suffixes", {
  tbl <- protocol_table_C0()

  labels <- purrr::pmap_chr(
    list(tbl$wavelength, tbl$stim_nd, tbl$stim_index, tbl$protocol_id),
    make_stimresp_label
  )

  # Block 1 (stim_index 1-7): 650A
  expect_true(all(grepl("^650A ", labels[1:7])))

  # Block 6 (stim_index 36-42): 650B
  expect_true(all(grepl("^650B ", labels[36:42])))

  # Other wavelengths have no A/B suffix
  non_650_idx <- which(tbl$wavelength != "650")
  expect_true(all(!grepl("^650[AB]", labels[non_650_idx])))
})

test_that("C0 sweep-to-stim mapping matches Origin for key positions", {
  # Validated against temp_file/comparison_results_export/protocol_compare.csv
  # Each sweep_in_protocol maps to exactly one (stim_index, stim_label)
  sweep_tbl <- expand_protocol_repeats(protocol_table_C0(), repeats_per_stim = 4L)

  sweep_labels <- purrr::pmap_chr(
    list(sweep_tbl$wavelength, sweep_tbl$stim_nd,
         sweep_tbl$stim_index, sweep_tbl$protocol_id),
    make_stimresp_label
  )

  # First 4 sweeps → stim_index 1 → "650A 4.0"
  expect_true(all(sweep_labels[1:4] == "650A 4.0"))
  expect_true(all(sweep_tbl$stim_index[1:4] == 1L))

  # Sweeps 141-144 → stim_index 36 → "650B 4.5"  (first rep of second 650 block)
  expect_true(all(sweep_labels[141:144] == "650B 4.5"))
  expect_true(all(sweep_tbl$stim_index[141:144] == 36L))

  # Last 4 sweeps (277-280) → stim_index 70 → "370 2.0"
  expect_true(all(sweep_labels[277:280] == "370 2.0"))
  expect_true(all(sweep_tbl$stim_index[277:280] == 70L))
})

test_that("C0 protocol table has required columns for downstream processing", {
  tbl <- protocol_table_C0()

  required_cols <- c(
    "protocol_id", "protocol_variant", "stim_type",
    "stim_index", "block_index", "within_block_index",
    "wavelength", "stim_nd"
  )

  expect_true(all(required_cols %in% names(tbl)))
  expect_equal(tbl$protocol_id[1], "C0")
  expect_equal(tbl$protocol_variant[1], "p4")
})

test_that("C0 expanded sweep table preserves all protocol metadata columns", {
  tbl  <- protocol_table_C0()
  swps <- expand_protocol_repeats(tbl, repeats_per_stim = 4L)

  # tech_rep runs 1-4 for every stim_index
  expect_equal(sort(unique(swps$tech_rep)), 1:4)

  # sweep_in_protocol is 1:280
  expect_equal(swps$sweep_in_protocol, seq_len(280L))

  # wavelength and stim_nd are preserved
  expect_true("wavelength" %in% names(swps))
  expect_true("stim_nd"    %in% names(swps))
})
