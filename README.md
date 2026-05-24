# Replication code for the RESTat streaming-curve paper

This repository contains the R code and saved outputs used for the simulations,
critical values, contaminated-training checks, and S&P 500 empirical application in
the paper on online monitoring of functional representations.

The repository is organized so that the main simulation study, supplementary
monitoring statistics, and empirical application can be inspected separately.

## Repository structure

### `main_statistics/`

R code for the main HAC, SSMS, and RSMS KS and weighted-CvM monitoring
statistics.

- `drivers/`: scripts for null simulations, power simulations, summaries, plots,
  and paper tables
- `R_modules/`: shared functions for data generation, FPCA, monitoring
  statistics, critical-value lookup, and utilities
- `critical_value_engine/`: Monte Carlo code for the main KS and weighted-CvM
  critical values

### `alternative_detectors/`

R code for the supplementary CUSUM and MOSUM monitoring statistics reported in
the appendix.

- `drivers/generate_alt_detector_critical_values.R`: driver for supplementary
  critical-value simulations
- `critical_value_engine/`: code for Page-CUSUM, weighted-CUSUM, and multiscale
  MOSUM critical values

### `contaminated_training/`

R code for the mild training-sample contamination experiment.

- `drivers/run_contaminated_training_streamingcurve_simulation.R`: runs the
  simulation
- `drivers/summarize_contaminated_training_streamingcurve.R`: summarizes the
  results
- `drivers/make_contaminated_training_plots.R`: creates the contamination figures
- `R_modules/`: shared functions for this simulation

### `empirical/`

Saved empirical results, figures, and scripts for the S&P 500 application.

- `data/`: saved empirical result tables
- `outputs/`: rebuilt empirical tables in CSV and LaTeX format
- `figures/`: empirical figures used in the paper
- `updated_code/`: scripts for rebuilding the empirical tables and summary
  figures
- `legacy_source/`: preserved source scripts for the original SPX preprocessing
  workflow

The raw 1-minute S&P 500 TickData file is licensed and is not included in this
public repository. To rerun the SPX figure script, place the licensed file at
`empirical/legacy_source/data/SPX.csv`.

### `notes/`

Short notes on the repository contents.

## Suggested commands

Run commands from the repository root unless noted otherwise.

To rebuild the main simulation tables and figures:

```bash
Rscript main_statistics/drivers/generate_paper_tables_figures_tex.R
```

To generate supplementary CUSUM and MOSUM critical values:

```bash
Rscript alternative_detectors/drivers/generate_alt_detector_critical_values.R
```

To rerun the contaminated-training simulation:

```bash
Rscript contaminated_training/drivers/run_contaminated_training_streamingcurve_simulation.R
```

To rebuild the empirical tables and available SPX figures:

```bash
Rscript empirical/updated_code/run_empirical_bundle_refresh.R
```

The empirical figure script requires the licensed `SPX.csv` file described above.
The saved empirical tables in `empirical/data/` and `empirical/outputs/` can be
inspected without that raw TickData file.

## Files not included

The repository does not include the licensed TickData input
`empirical/legacy_source/data/SPX.csv`. It also omits archived duplicates,
temporary working copies, and unrelated older project versions.
