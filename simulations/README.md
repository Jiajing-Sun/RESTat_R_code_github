# Simulation Code

This folder contains current R code for the Monte Carlo studies in the paper and appendices.

## Subfolders

- `main/`
  Main null, power, size-adjusted power, average detection delay, HAC sensitivity, and multiscale tuning scripts.
- `contaminated_training/`
  Simulations for mild training contamination and training-window robustness designs.
- `diagnostics/`
  Diagnostic scripts, including the RSMS off-diagonal diagnostic.

## Data and outputs

No raw simulation output, summary CSV, cached critical values, or generated figures are included. Running the scripts
will create local output folders.

## Typical entry points

From `simulations/main/`:

```bash
Rscript run_null_streamingcurve_simulation.R
Rscript run_power_streamingcurve_simulation.R
Rscript summarize_streamingcurve_simulations.R
Rscript build_latest_simulation_section.R
```

From `simulations/contaminated_training/`:

```bash
Rscript run_contaminated_training_streamingcurve_simulation.R
Rscript summarize_contaminated_training_streamingcurve.R
Rscript build_contaminated_training_paper_assets.R
```
