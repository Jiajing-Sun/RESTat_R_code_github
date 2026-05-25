# Empirical Analysis Code

This folder contains the current empirical-analysis scripts for the S&P 500 intraday return-curve application.

Only the current empirical scripts are included.

## Data

The TickData 1-minute S&P 500 file `SPX.csv` is not included. To run the empirical scripts, place the licensed data
file at:

```text
empirical_analysis/data/SPX.csv
```

or pass the path explicitly:

```bash
Rscript run_empirical_bundle_refresh.R --data_file=/path/to/SPX.csv
```

Individual scripts also read the `SPX_DATA_FILE` environment variable.

## Current scripts

- `prepare_spx_summary_figures.R`
- `build_intraday_curve_regime_figure.R`
- `build_intraday_surface_detection_figure.R`
- `rebuild_empirical_outputs_from_saved_results.R`
- `compute_decision_window_diagnostics.R`
- `compute_q_robustness_diagnostics.R`
- `compute_training_window_robustness.R`
- `run_empirical_bundle_refresh.R`

Generated figures, CSV outputs, and saved empirical result tables are not included.

`rebuild_empirical_outputs_from_saved_results.R` additionally uses optional cached empirical-result files
`data/empirical_main_results.csv` and `data/empirical_alt_results.csv`. The public bundle does not include those files;
`run_empirical_bundle_refresh.R` skips that cache-dependent step when they are absent.
