# Empirical bundle for the current paper version

This folder collects the empirical material that matches the current paper version:

- `Online_Monitoring_via_Streaming_Curves-Final.tex`
- `Online_Monitoring_of_Structural_Change_in_Functional_Time_Series_via_Self_Normalization_Online_Supplement.tex`

## Folder layout

### `legacy_source/`
Recovered legacy empirical material from the uploaded `empircal/` folder.

- `legacy_source/data/README.md`
  Documents the licensed TickData input `SPX.csv` expected by the empirical scripts. The raw TickData file is not redistributed in this public bundle.
- `legacy_source/code/`
  Legacy SPX scripts and supporting functions, including:
  - `sp500-test-20260101.R`
  - `sp500-test.R`
  - `Method.R`
  - `Method-fct20260102.R`
  - `genData.R`
  - `spx-plot.R`
  - `spx-plot_RR.R`

These scripts are preserved because they are part of the empirical project history. They are useful for understanding the older workflow and for recovering the original SPX preprocessing logic.

### `data/`
Saved empirical result data used by the current paper version.

- `empirical_main_results.csv`
- `empirical_alt_results.csv`
- `empirical_main_detection_rank.csv`
- `empirical_alt_detection_rank.csv`
- `empirical_narrative_notes.txt`

These are the preserved current-paper empirical outputs extracted from the project workspace.

### `outputs/`
Paper-facing empirical tables for the current manuscript.

This folder contains:

- rebuilt main KS and weighted-CvM empirical tables
- a full alternative-detector table copied from the historical outputs
- a retained alternative-detector table matching the current supplement emphasis

### `figures/`
Current empirical figures used by the paper:

- `SPX_2020_Min_Close.png`
- `SPX_2020_DailyLogReturn.png`

### `updated_code/`
Current helper scripts for the revised paper version.

- `prepare_spx_summary_figures.R`
  Regenerates the two SPX summary figures when the licensed `legacy_source/data/SPX.csv` file is available locally
- `rebuild_empirical_outputs_from_saved_results.R`
  Rebuilds the current empirical tables, rankings, and notes from the preserved saved results in `data/`
- `run_empirical_bundle_refresh.R`
  Convenience entry point that runs both steps above

## Important note

The recovered legacy SPX scripts do not by themselves constitute a complete end-to-end driver for the current revised empirical section with the expanded KS / weighted-CvM / retained alternative-detector reporting structure. For that reason, this bundle keeps both:

- the recovered legacy source base
- the current paper-matched rebuild layer in `updated_code/`

This is intentional and reflects the actual project history more faithfully than pretending that one legacy script alone fully reproduces the current manuscript.

The raw 1-minute S&P 500 input `SPX.csv` comes from TickData and is not included in the public repository. To rerun the SPX figure script, place the licensed file at `legacy_source/data/SPX.csv`.

## Suggested usage

From this folder, run:

```bash
Rscript updated_code/run_empirical_bundle_refresh.R
```

This will:

1. regenerate the current SPX figures in `figures/`
2. rebuild the current empirical tables and rankings in `outputs/`
