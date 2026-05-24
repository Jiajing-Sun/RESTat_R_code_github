# Streaming-curve critical values (fresh simulation)

This folder generates critical values **from scratch** for the streaming-curve paper.
It does **not** warm-start from older tables and does **not** top up missing entries.

## Statistics covered

- Shao's KS (`SSMS`, `KS`)
- Adjusted-range based KS (`RSMS`, `KS`)
- HAC-based KS (`HAC`, `KS`)
- Shao's weighted CvM (`SSMS`, `CvM`)
- Adjusted-range based weighted CvM (`RSMS`, `CvM`)
- HAC-based weighted CvM (`HAC`, `CvM`)

## Default simulation design

- `q = 1:30`
- `T = 1, 2, 5, 10, Inf`
- `gamma = c(0, 0.15)` for KS statistics
- CvM weights: `U`, `Early`, `Mid`, `Late`
- significance levels: `0.10`, `0.05`, `0.01`
- finite-horizon replications: `10000`
- open-end replications: `5000`

## How to run

### RStudio / interactive R

Open `run_generate_streamingcurve_cv.R` and run it.

### Rscript / terminal

```r
Rscript run_generate_streamingcurve_cv.R
```

You may override the project root or the number of workers from the command line:

```r
Rscript run_generate_streamingcurve_cv.R --root=/path/to/streaming_curve_cv_codes --ncores=8
```

On Windows Server:

```r
Rscript.exe run_generate_streamingcurve_cv.R --ncores=8
```

## Output files

The script writes:

- `outputs/critical_values_all.csv`
- `outputs/critical_values_all_weights.csv`
- `outputs/run_summary.txt`

If files with the same names already exist, they are backed up with a timestamp and then overwritten.

## Notes

- The script is designed to be cross-platform and does not require manual `setwd()`.
- Parallel execution uses a PSOCK cluster so it works on Windows.
- The `critical_values.R` helper provides loading and lookup functions for the generated CSV files.
