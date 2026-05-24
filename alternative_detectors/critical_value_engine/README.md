
Alternative-detector critical-value simulation (fast fixed version)
=================================================================

Why the previous script appeared to hang
---------------------------------------
The earlier version scanned Page-CUSUM / weighted-CUSUM over very large endpoint-lag grids,
repeating that scan separately for every q=1,...,30 and for every standardizer. That is not a
literal deadlock, but it is computationally explosive in plain R and can easily run for a day
without visible progress.

What changed in this version
----------------------------
1. Batched progress reporting: you now see batch-by-batch progress.
2. Deterministic scan thinning: the Brownian path is still generated on the requested grids,
   but benchmark detector scans are evaluated on a thinned endpoint grid by default.
3. Much faster q-handling: all q=1,...,30 are computed jointly instead of rerunning the same
   scan 30 times.
4. Fresh-run output behavior remains the same.

Important note
--------------
This version is designed to be feasible in plain R. For the benchmark detector families
(Page-CUSUM, weighted-CUSUM, MOSUM, multiscale MOSUM), the default settings therefore use
thinned scan grids rather than exact all-endpoint/all-lag scans.

Main script
-----------
run_generate_alt_detector_cv.R

Recommended first run on Windows Server
--------------------------------------
Rscript run_generate_alt_detector_cv.R --ncores=12

Outputs
-------
outputs/critical_values_alt_detectors.csv
outputs/run_summary.txt
