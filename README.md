
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ZEUS

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

ZEUS is an R package for analyzing electroretinography (ERG) recordings
stored in Axon Binary Format (`.abf`) files. It is designed around the
ZEUS stimulus protocols used in the ZENV workflow and provides a
reproducible pipeline for importing raw recordings, mapping sweeps to
protocol metadata, generating averaged stimulus-response traces,
extracting peak summaries, validating outputs against reference exports,
and exporting analysis-ready tables.

The package supports both:

- `C0`: spectral recordings split across wavelength blocks and
  neutral-density (ND) levels
- `C1`: white-light recordings organized as repeated ND sweeps

## Why use ZEUS

ZEUS turns a raw ABF recording into a structured, reproducible analysis
object. Instead of re-labeling sweeps manually or relying on
spreadsheet-only post-processing, you can keep the full workflow inside
R and produce outputs that are easier to validate, version, and reuse.

Core capabilities include:

- ABF import with `readABF`
- sweep-to-protocol labeling for ZEUS `C0` and `C1`
- replicate-aware averaging into 70-condition stimulus-response traces
- baseline correction, smoothing, and noisy-trace filtering
- waveform plotting for mean, spectral, and intensity-response summaries
- A-wave, B-wave, and D-wave peak summary tables
- CSV bundle and Excel workbook export helpers
- validation utilities for comparison against Origin StimResp exports
- a Shiny app for interactive review and export

## Installation

ZEUS is currently a development package and can be installed from
GitHub.

``` r
# install.packages("pak")
pak::pak("lo2252/zeus")
```

You can also install it with `remotes`:

``` r
# install.packages("remotes")
remotes::install_github("lo2252/zeus")
```

ZEUS depends on R 4.1 or later.

## Package dependencies

The authoritative dependency list lives in `DESCRIPTION`, but for
convenience the main package dependencies are listed here as well.

Required runtime packages:

- `bslib`
- `dplyr`
- `DT`
- `ggplot2`
- `pkgload`
- `purrr`
- `readABF`
- `rlang`
- `scales`
- `shiny`
- `signal`
- `stringr`
- `tibble`
- `tidyr`
- `zoo`

Suggested packages used for development or testing:

- `rstudioapi`
- `testthat`

If you install ZEUS from GitHub with `pak::pak()` or
`remotes::install_github()`, these package dependencies are installed
automatically.

## Included example data

The package ships with example ABF files in `inst/extdata` so you can
try the workflow without providing your own recordings immediately.

| File           | Intended protocol | Description                       |
|----------------|-------------------|-----------------------------------|
| `26225004.abf` | `C1`              | White-light ND sweep example      |
| `26225005.abf` | `C0`              | Spectral wavelength-block example |

Both example recordings contain 280 sweeps and are used in package tests
and validation scripts.

## Quick start

The main entry point is `zeus_read_abf()`, which reads an ABF file, maps
sweeps to the requested ZEUS protocol, and returns a structured
`zeus_stimresp` object.

``` r
library(ZEUS)

abf_path <- system.file("extdata", "26225004.abf", package = "ZEUS")

x <- zeus_read_abf(
  path = abf_path,
  protocol = "C1"
)

class(x)
names(x)
```

The returned object typically includes:

- `raw`: the imported ABF object
- `traces_280`: sweep-level ERG traces labeled with protocol metadata
- `traces_70`: averaged stimulus-response traces collapsed to 70
  conditions
- `photocell`: photocell traces when present
- `stimresp_qc`: quality-control summaries
- `stimresp_settings`: the processing settings used to build the object

## Protocol overview

ZEUS includes protocol builders for both supported acquisition modes:

``` r
protocol_table_C0()
protocol_table_C1()
```

At a high level:

- `C0` defines 70 spectral conditions arranged as 10 wavelength blocks x
  7 ND levels
- `C1` defines 70 white-light conditions arranged as 10 runs x 7 ND
  levels
- `expand_protocol_repeats()` expands those 70 conditions into the
  280-sweep order used in the raw recordings

This protocol metadata is joined directly onto sweep-level waveform
tables so downstream summaries remain traceable back to the acquisition
structure.

## Typical analysis workflow

### 1. Import and label a recording

``` r
x <- zeus_read_abf(
  path = system.file("extdata", "26225005.abf", package = "ZEUS"),
  protocol = "C0",
  align_to_stimulus = "protocol"
)
```

### 2. Inspect the averaged traces

``` r
head(x$traces_70)
unique(x$traces_70$stim_label)
```

### 3. Plot mean waveforms

``` r
zeus_plot_mean_waveform(x)
```

For spectral recordings, you can generate a faceted wavelength-block
view:

``` r
zeus_plot_spectral_waveform(x)
```

For intensity-response summaries across ND levels:

``` r
zeus_plot_intensity_response(x)
```

### 4. Summarize peak statistics

``` r
peak_summary <- zeus_summarize_peak_statistics(x)

peak_summary$key_statistics
peak_summary$by_nd
peak_summary$by_wavelength
```

The summary output is designed for reporting and export. It includes
combined tables for A-wave, B-wave, and D-wave measurements, plus
ND-level and wavelength-level aggregations where applicable.

## Exporting results

ZEUS can write the major outputs from a processed object to either a set
of CSV files or a single Excel workbook.

### Export a CSV bundle

``` r
zeus_export_csv_bundle(
  x,
  csv_path = "exports/example_run"
)
```

This writes component files such as:

- `example_run_raw.csv`
- `example_run_traces_280.csv`
- `example_run_traces_70.csv`
- `example_run_stimresp_qc.csv`
- `example_run_stimresp_settings.csv`
- `example_run_peak_statistics.csv`

### Export an Excel workbook

``` r
zeus_export_excel_workbook(
  x,
  xlsx_path = "exports/example_run.xlsx"
)
```

## Validation against reference exports

One of the package goals is reproducibility against established
ZEUS/Origin workflows. The repository includes validation scripts in
`validation/` for both example protocols:

- `validation/validate_c0.R`
- `validation/validate_c1.R`

These scripts:

- load the package from source
- import the example ABF recording with `zeus_read_abf()`
- compare protocol labels and sweep mappings against Origin reference
  exports
- compare waveform outputs point-by-point
- compute agreement summaries such as correlation, mean absolute
  difference, and RMSE

For custom validation pipelines, ZEUS also exports helper functions such
as:

``` r
zeus_compare_waveforms(...)
zeus_summarize_protocol_validation(...)
zeus_validate_protocol_agreement(...)
zeus_validate_response_agreement(...)
```

## Shiny application

ZEUS includes an interactive Shiny application for users who want a
guided workflow for import, visualization, statistics review, and
export.

``` r
run_zeus_app()
```

The app is useful when you want to:

- review traces before exporting results
- compare waveform views interactively
- inspect summary tables without writing additional R code
- prepare export bundles for downstream analysis

## Development status

ZEUS is under active development and the current lifecycle is marked
experimental. Public interfaces may continue to evolve as the package
grows, but the core import, protocol labeling, plotting, validation, and
export workflow is already in place.

## Contributing

Issues and pull requests are welcome. If you are extending the package,
edit `README.Rmd` rather than `README.md`, then regenerate the Markdown
README so the rendered file stays in sync.
