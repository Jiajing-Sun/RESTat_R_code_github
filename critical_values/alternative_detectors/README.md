# Alternative-Detector Critical Values

This folder contains current code for simulating critical values for alternative detector families.

- `run_generate_alt_detector_cv.R` is the direct alternative-detector critical-value engine.
- `generate_alt_detector_critical_values.R` is the current simulation-bundle driver for alternative-detector critical values.
- `R/` contains both the alternative-detector critical-value engine helpers and the current support modules needed by the driver.

Generated critical-value tables are not included. The default command writes finite-horizon tables to `outputs/`:

```bash
Rscript run_generate_alt_detector_cv.R --ncores=4
```

The public default excludes `T=Inf`; add `--include_openend=true` only for separate open-end experiments.
