## Test environments

- Local macOS (aarch64-apple-darwin20), R 4.5.3

## R CMD check results

- `devtools::test()` passed with 130 tests, 0 failures, 0 warnings, 0 skips.
- `R CMD check --as-cran ZEUS_0.1.1.tar.gz` completed with 0 errors, 0
  warnings, and 1 NOTE.

## Notes

- This is a new submission.
- Local `R CMD check --as-cran ZEUS_0.1.1.tar.gz` completed with the expected
  NOTE for a new submission.
- The package includes bundled example ABF files in `inst/extdata` for testing,
  validation, and demonstration workflows.
