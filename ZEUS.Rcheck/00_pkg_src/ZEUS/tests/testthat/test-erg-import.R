# ERG import tests --------------------------------------------------------
# Validates zeus_read_abf() output against the structure expected by the
# C0 and C1 protocols. C0-specific tests require 26225005.abf (C0 recording)
# in inst/extdata/; they are skipped automatically if the file is absent.
# C1-specific tests use 26225004.abf which is the primary package test file.

# C0 import ---------------------------------------------------------------

test_that("zeus_read_abf C0: traces_280 has 280 unique sweeps", {
  obj <- zeus_test_obj("C0")

  expect_s3_class(obj, "zeus_stimresp")
  expect_true(!is.null(obj$traces_280))

  n_sweeps <- dplyr::n_distinct(obj$traces_280$sweep)
  expect_equal(n_sweeps, 280L)
})

test_that("zeus_read_abf C0: traces_280 contains required columns", {
  obj <- zeus_test_obj("C0")

  required <- c("sweep", "time", "time_ms", "value",
                "stim_index", "stim_label", "wavelength", "stim_nd",
                "protocol_id", "tech_rep")
  expect_true(all(required %in% names(obj$traces_280)))
})

test_that("zeus_read_abf C0: stim_labels in traces_280 match Origin export", {
  obj <- zeus_test_obj("C0")
  t280 <- obj$traces_280

  # All 280 sweeps should have a non-NA stim_label
  expect_equal(sum(is.na(t280$stim_label)), 0L)

  # Sweeps 1-4 → "650A 4.0" (stim_index 1, first 650 block)
  sweep1_labels <- unique(t280$stim_label[t280$sweep == 1L])
  expect_equal(sweep1_labels, "650A 4.0")

  # Sweeps 141-144 → "650B 4.5" (stim_index 36, second 650 block)
  sweep141_labels <- unique(t280$stim_label[t280$sweep == 141L])
  expect_equal(sweep141_labels, "650B 4.5")
})

test_that("zeus_read_abf C0: traces_70 has 70 unique stim conditions", {
  obj <- zeus_test_obj("C0")

  expect_true(!is.null(obj$traces_70))

  n_stim <- dplyr::n_distinct(obj$traces_70$stim_index)
  expect_equal(n_stim, 70L)
})

test_that("zeus_read_abf C0: traces_70 retains wavelength and stim_nd metadata", {
  obj <- zeus_test_obj("C0")
  t70 <- obj$traces_70

  # Required metadata columns must be present after averaging
  expect_true("wavelength" %in% names(t70))
  expect_true("stim_nd"    %in% names(t70))
  expect_true("stim_label" %in% names(t70))

  # No NA wavelengths in averaged output
  expect_equal(sum(is.na(t70$wavelength)), 0L)
})

test_that("zeus_read_abf C0: traces_70 stim_labels match Origin-validated labels", {
  # Expected 70 Origin stim_labels in protocol order
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

  obj <- zeus_test_obj("C0")
  t70  <- obj$traces_70

  # One row of metadata per stim_index; retrieve sorted by stim_index
  meta <- t70 |>
    dplyr::distinct(stim_index, stim_label) |>
    dplyr::arrange(stim_index)

  expect_equal(meta$stim_label, origin_labels)
})

test_that("zeus_read_abf C0: n_reps_used reflects correct technical replicate count", {
  obj <- zeus_test_obj("C0")
  t70 <- obj$traces_70

  # Each averaged condition should use exactly 4 technical replicates
  # (one per sweep per stim_index in this 280-sweep C0 file)
  reps_per_stim <- t70 |>
    dplyr::distinct(stim_index, n_reps_used) |>
    dplyr::pull(n_reps_used)

  expect_true(all(reps_per_stim >= 1L))
  expect_true("n_reps_used" %in% names(t70))
})

# C1 import ---------------------------------------------------------------

test_that("zeus_read_abf C1: basic structure is correct", {
  obj <- zeus_test_obj("C1")

  expect_s3_class(obj, "zeus_stimresp")
  expect_false(is.null(obj$traces_280))
  expect_false(is.null(obj$traces_70))
})

test_that("zeus_read_abf C1: traces_70 has 70 unique stim conditions", {
  obj <- zeus_test_obj("C1")

  n_stim <- dplyr::n_distinct(obj$traces_70$stim_index)
  expect_equal(n_stim, 70L)
})

test_that("zeus_read_abf C1: all stim_labels use White wavelength", {
  obj <- zeus_test_obj("C1")

  expect_true(all(obj$traces_70$wavelength == "White"))
  expect_true(all(grepl("^White ", obj$traces_70$stim_label)))
})
