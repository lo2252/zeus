# validation/

Scripts for validating ZEUS C0 and C1 protocol transformations against Origin
StimResp exports.

## Scripts

| File              | Protocol | ABF source                      | Origin reference                              |
|-------------------|----------|---------------------------------|-----------------------------------------------|
| `validate_c0.R`   | C0       | `temp_file/26225005.abf`         | `temp_file/26225005_origin_export_with_d_wave.xlsx` |
| `validate_c1.R`   | C1       | `inst/extdata/26225004.abf`      | `temp_file/26225004_origin_export_with_d_wave.xlsx` |

## Usage

Run from the package root directory (where `DESCRIPTION` lives):

```r
# In an R session:
setwd("<path-to-zeus>")
source("validation/validate_c0.R")
source("validation/validate_c1.R")
```

Or from the command line:
```bash
Rscript validation/validate_c0.R
Rscript validation/validate_c1.R
```

Each script:
1. Loads the ZEUS package via `pkgload::load_all()`.
2. Imports the ABF file with `zeus_read_abf()`.
3. Imports the Origin StimResp sheet from the xlsx file.
4. Compares the 70-condition protocol labels (by stim_index) to Origin.
5. Compares the 280-sweep mapping to Origin.
6. Computes mean ERG traces from ZEUS and compares point-by-point to Origin
   StimResp traces (RMSE, mean |diff|, Pearson correlation).
7. Prints a summary table and stores results in `results_c0` / `results_c1`.

## Interpreting results

- **Protocol label matches 70/70** → the condition ordering is identical to Origin.
- **Sweep label matches 280/280** → the sweep-to-condition mapping is correct.
- **Trace correlation ≈ 1.0** → ZEUS mean traces are numerically identical to
  Origin StimResp traces (after baseline correction and smoothing; small
  differences are expected due to processing differences).
