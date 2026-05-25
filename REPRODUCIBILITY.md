# Reproducibility Guide

This repository separates quick code checks from full Monte Carlo reproduction.

## Lightweight Checks

Run:

```sh
Rscript requirements.R --check
./run_smoke_tests.sh
```

The smoke tests parse all R files, check for private artifacts, and run a tiny finite-horizon critical-value job in a
temporary output directory. They are intended to catch path, syntax, and packaging problems after a fresh clone.

## Critical Values

The full critical-value settings are:

- main finite-horizon critical values: `nrep_finite = 10000`
- open-end critical values, if explicitly requested: `nrep_openend = 5000`
- public default horizons: `T = 1, 2, 5, 10`

Generate the public finite-horizon tables with:

```sh
Rscript critical_values/main/run_generate_streamingcurve_cv.R --ncores=4
Rscript critical_values/alternative_detectors/run_generate_alt_detector_cv.R --ncores=4
```

The resulting CSV files are written under `critical_values/*/outputs/` and are ignored by git.

## Simulations

After critical values are generated, run the Monte Carlo drivers under `simulations/main/` and
`simulations/contaminated_training/`. The full simulation setting is `nsim = 1000` per design unless a script
explicitly states otherwise.

## Empirical Analysis

The empirical scripts require licensed TickData input. See `DATA_ACCESS.md`. Some scripts also rebuild tables from
saved empirical-result caches (`empirical_main_results.csv` and `empirical_alt_results.csv`) when those files are available.
The public bundle does not include those caches; the empirical refresh driver skips cache-dependent steps when they are
missing.
