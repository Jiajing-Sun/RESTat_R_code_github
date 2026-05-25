# Critical-Value Simulation Code

This folder contains current R code for critical-value simulation. Public defaults are finite-horizon only.

## Subfolders

- `main/`
  Critical-value engine for the main KS and weighted-CvM monitor families.
- `alternative_detectors/`
  Critical-value code for alternative detector shapes such as Page-CUSUM, weighted CUSUM, and multiscale MOSUM.

## Data and outputs

No generated critical-value CSV files are included. Running the scripts will create output files locally under
`main/outputs/` or `alternative_detectors/outputs/`.

## Typical entry points

From `critical_values/main/`:

```bash
Rscript run_generate_streamingcurve_cv.R --ncores=4
```

From `critical_values/alternative_detectors/`:

```bash
Rscript run_generate_alt_detector_cv.R --ncores=4
```

or, for the simulation-support generator used in the current code bundle:

```bash
Rscript generate_alt_detector_critical_values.R --ncores=4
```

The scripts default to `T = 1, 2, 5, 10`. Open-end helper routines are retained for experiments but are not run unless
`--include_openend=true` is supplied explicitly.
