# Paper-Matched R Code Bundle

This folder collects the R source files that match the current paper version:

- `Online_Monitoring_via_Streaming_Curves-Final.tex`
- `Online_Monitoring_of_Structural_Change_in_Functional_Time_Series_via_Self_Normalization_Online_Supplement.tex`

The goal is to keep one clean code bundle for the current main paper and online supplement, without mixing in older duplicates, archived copies, or temporary variants such as files with names ending in ` 2.R`.

## Folder structure

### `main_statistics/`
Code for the main KS and weighted-CvM procedures studied in the main paper.

- `drivers/`
  - `run_null_streamingcurve_simulation.R`: null-size simulations for the main monitoring family
  - `run_power_streamingcurve_simulation.R`: power simulations for the main monitoring family
  - `summarize_streamingcurve_simulations.R`: builds null / power / size-adjusted power / ADL summaries
  - `make_streamingcurve_plots.R`: produces the simulation figures
  - `generate_paper_tables_figures_tex.R`: project-level entry point for regenerating the simulation tables / figures / LaTeX
  - `build_latest_simulation_section.R`: writes the simulation material used in the current paper version
- `R_modules/`
  Shared modules used by the main simulation pipeline, including FPCA, data generation, detector definitions, method catalogs, and utility functions
- `critical_value_engine/`
  Monte Carlo code for the main-statistics critical values used in the paper

This block corresponds to:

- the main simulation section in `Online_Monitoring_via_Streaming_Curves-Final.tex`
- the main-statistics critical-value appendix in `Online_Monitoring_via_Streaming_Curves-Final.tex`

### `alternative_detectors/`
Code for the retained alternative detector shapes that appear in the online supplement.

- `drivers/`
  - `generate_alt_detector_critical_values.R`: project-facing driver for the retained alternative detector critical values
- `critical_value_engine/`
  Low-level code for Page-CUSUM, weighted-CUSUM, and multiscale-MOSUM critical-value generation

This block corresponds to:

- the alternative-detector critical-value appendix in the online supplement
- the retained alternative-detector simulation appendix in the online supplement

Note: the supplement text keeps Page-CUSUM as a secondary comparison and emphasizes weighted-CUSUM and multiscale-MOSUM. Full-CUSUM and single-bandwidth MOSUM were computed in the project history but are not the retained focus of the current supplement narrative.

### `contaminated_training/`
Code for the mild Phase-I contamination study.

- `drivers/`
  - `run_contaminated_training_streamingcurve_simulation.R`: runs the contamination experiment
  - `summarize_contaminated_training_streamingcurve.R`: summarizes the contamination results
  - `make_contaminated_training_plots.R`: generates contamination-study figures
  - `build_contaminated_training_paper_assets.R`: writes paper-facing tables / figures / LaTeX-ready assets
- `R_modules/`
  Shared modules for contamination design, detectors, FPCA, utilities, and simulation support

This block corresponds to:

- the contamination-robustness subsection in `Online_Monitoring_via_Streaming_Curves-Final.tex`
- the contamination appendix material in the online supplement

### `empirical/`
This subfolder now contains both:

- the recovered legacy empirical source and data from `empircal/`
- an updated paper-matched rebuild layer for the current empirical tables and figures

The empirical bundle is organized as follows:

- `legacy_source/data/`
  - placeholder documentation for the licensed TickData input `SPX.csv`, which is not redistributed in this public bundle
- `legacy_source/code/`
  - recovered legacy empirical scripts such as `sp500-test-20260101.R`, `sp500-test.R`, `Method.R`, and the original SPX plotting files
- `data/`
  - saved empirical result tables used by the current paper version
- `outputs/`
  - paper-facing empirical tables rebuilt from the saved results
- `figures/`
  - empirical figures used by the current paper version
- `updated_code/`
  - current helper scripts that regenerate the SPX summary figures and rebuild the empirical tables / rankings from the saved results

Important distinction:

The recovered legacy SPX scripts predate the final expanded method family in the revised paper and therefore should be treated as the historical empirical source base, not as a complete end-to-end driver for the current `Final.tex` and supplement. The new `updated_code/` layer bridges that gap by rebuilding the current paper-facing empirical assets from the preserved results. Regenerating the SPX summary figures requires the licensed TickData input `SPX.csv`, which should be placed locally under `empirical/legacy_source/data/`.

See `empirical/README.md` for the exact layout and usage.

### `notes/`
Short documentation files for provenance and bundle status.

## What was intentionally excluded

This bundle excludes:

- archived duplicates under nested historical folders
- temporary copies such as `* 2.R`
- diagnostic or audit scripts that are not needed to reproduce the current paper-facing results
- unrelated project versions such as `Online_Monitoring_via_Streaming_Curves_v1/`
- the licensed TickData raw input `empirical/legacy_source/data/SPX.csv`

## Original source locations

The files in this bundle were copied from the current project tree, mainly from:

- `fresh_streaming_curve_sim_codes_副本/`
- `streaming_curve_contam_training_sim_codes/`
- the project root `build_latest_simulation_section.R`

They were selected to match the current `Final.tex` and online supplement, rather than older manuscript versions.

## Suggested usage

If you want to regenerate the current paper-facing simulation assets, start from:

```bash
Rscript main_statistics/drivers/generate_paper_tables_figures_tex.R
```

If you want to regenerate the retained alternative-detector critical values, start from:

```bash
Rscript alternative_detectors/drivers/generate_alt_detector_critical_values.R
```

If you want to rerun the contaminated-training study, start from:

```bash
Rscript contaminated_training/drivers/run_contaminated_training_streamingcurve_simulation.R
```

## Matching status

This bundle is intended to match the current paper structure:

- main paper: HAC / SSMS / RSMS KS and weighted-CvM, plus the retained contamination study
- online supplement: retained alternative detector shapes, their critical values, additional simulation evidence, contamination robustness, and supplementary empirical outputs

For the empirical block, the bundle now contains:

- recovered legacy raw data and scripts
- saved current-paper empirical outputs
- updated rebuild scripts for the current paper-facing empirical figures and tables
