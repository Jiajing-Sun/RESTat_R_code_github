# Main Simulation Code

This folder contains the current main simulation drivers and helper modules.

- Driver scripts are in this folder.
- `R/` contains the current modules used by the drivers.

The paper reports 1000 Monte Carlo replications for the simulation studies. Critical-value simulations use the
replication counts specified in the corresponding critical-value scripts.

Generated `outputs/`, `figures/`, and local critical-value CSV files are not included. Before running the full
simulation drivers, generate the needed critical values under `../../critical_values/*/outputs/` or provide equivalent
CSV files in one of the lookup paths reported by `R/critical_values_lookup.R`.
