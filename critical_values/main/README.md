# Main Critical Values

This folder contains the current critical-value engine for the main KS and weighted-CvM monitors.

- `run_generate_streamingcurve_cv.R` is the entry script.
- `R/` contains helper functions for Brownian simulation, weight definitions, I/O, and path handling.

Generated critical-value tables are not included. The default command writes finite-horizon tables to `outputs/`:

```bash
Rscript run_generate_streamingcurve_cv.R --ncores=4
```

Use `--q_grid=1:5`, `--T_grid=1,2`, or `--nrep_finite=100` for lightweight local tests. The public default excludes
`T=Inf`; add `--include_openend=true` only for separate open-end experiments.
