# Empirical S&P 500 application

This folder contains the saved results, figures, and R scripts for the S&P 500
application in the paper.

## Folder layout

### `data/`

Saved empirical result files used to build the paper tables:

- `empirical_main_results.csv`
- `empirical_alt_results.csv`
- `empirical_main_detection_rank.csv`
- `empirical_alt_detection_rank.csv`
- `empirical_narrative_notes.txt`

### `outputs/`

Empirical tables rebuilt from the saved result files. The folder contains CSV and
LaTeX versions of the main KS, weighted-CvM, and supplementary monitoring
statistic tables.

### `figures/`

Empirical figures used in the paper:

- `SPX_2020_Min_Close.png`
- `SPX_2020_DailyLogReturn.png`

### `updated_code/`

Scripts for rebuilding the empirical outputs:

- `rebuild_empirical_outputs_from_saved_results.R`: rebuilds empirical tables,
  rankings, and notes from the saved results in `data/`
- `prepare_spx_summary_figures.R`: regenerates the two SPX summary figures when
  the licensed TickData file is available locally
- `run_empirical_bundle_refresh.R`: runs both steps above

### `legacy_source/`

Preserved SPX preprocessing scripts and supporting functions. These scripts are
included for transparency about the empirical data preparation. The standard
entry point for rebuilding the current empirical tables is
`updated_code/run_empirical_bundle_refresh.R`.

The raw 1-minute S&P 500 input `SPX.csv` comes from TickData and is not
redistributed in this public repository. To rerun the SPX figure script, place
the licensed file at:

```text
empirical/legacy_source/data/SPX.csv
```

## Suggested command

From the repository root, run:

```bash
Rscript empirical/updated_code/run_empirical_bundle_refresh.R
```

The saved empirical tables can be rebuilt without the raw TickData file. The SPX
summary figures require `SPX.csv`.
