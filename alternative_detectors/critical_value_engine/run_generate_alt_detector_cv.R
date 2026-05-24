
# ==============================================================
# run_generate_alt_detector_cv.R
# ============================================================== 
# Fresh simulation driver for alternative detector shapes after FPCA compression.
# Fast/plain-R version with thinned scan grids and batch progress.

bootstrap_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p1)) return(normalizePath(p1, winslash = "/", mustWork = FALSE))
    p2 <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  NULL
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
source_bootstrap_file <- function(rel_path) {
  candidates <- unique(c(file.path(dirname(bootstrap_script_path() %||% getwd()), rel_path), file.path(getwd(), rel_path)))
  for (p in candidates) if (file.exists(p)) { source(p, local = FALSE); return(invisible(p)) }
  stop("Could not source bootstrap helper: ", rel_path)
}
source_bootstrap_file(file.path("R", "project_paths.R"))
BOOTSTRAP_PATH <- bootstrap_script_path(); BOOTSTRAP_DIR <- if (!is.null(BOOTSTRAP_PATH)) dirname(BOOTSTRAP_PATH) else getwd(); ROOT <- resolve_project_root(default_start = BOOTSTRAP_DIR)
source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "alt_detector_weights.R"), local = FALSE)
source(file.path(ROOT, "R", "finite_alt_detector_critical_values.R"), local = FALSE)
source(file.path(ROOT, "R", "openend_alt_detector_critical_values.R"), local = FALSE)
source(file.path(ROOT, "R", "critical_values_io_alt.R"), local = FALSE)

# User controls
q_grid <- 1:30
T_grid <- c(1, 2, 5, 10, Inf)
gamma_vec <- c(0, 0.15)
alpha_levels <- c(0.10, 0.05, 0.01)
mosum_h_vec <- c(0.10, 0.20)
weighted_cusum_omega_names <- c("InvSqrt")
multiscale_h_sets <- list(H050_100_200 = c(0.05, 0.10, 0.20))
multiscale_scale_names <- c("Equal")
page_length_grid_size <- 40L
weighted_length_grid_size <- 40L
open_pair_lag_grid_size <- 40L
finite_eval_grid_size <- 250L
open_eval_grid_size <- 250L
exact_page_scan <- FALSE
exact_weighted_scan <- FALSE
nrep_finite <- 10000L
nrep_openend <- 5000L
grid_monitor <- 10000L
n_train_grid_open <- 1500L
n_open_grid_open <- 2000L
batch_size <- 100L
ridge <- 1e-10
range_floor <- 1e-8
args <- parse_named_args(); ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg); ncores <- core_plan$used
seed_finite <- 24681357L; seed_openend <- 97531L; verbose <- TRUE; backup_existing_outputs <- TRUE
OUTPUT_DIR <- file.path(ROOT, "outputs"); ensure_dir(OUTPUT_DIR)
summary_lines <- report_simulation_design_alt(root = ROOT, output_dir = OUTPUT_DIR, q_grid = q_grid, T_grid = T_grid,
  gamma_vec = gamma_vec, alpha_levels = alpha_levels, nrep_finite = nrep_finite, nrep_openend = nrep_openend,
  grid_monitor = grid_monitor, n_train_grid_open = n_train_grid_open, n_open_grid_open = n_open_grid_open,
  core_plan = core_plan, mosum_h_vec = mosum_h_vec, weighted_cusum_omega_names = weighted_cusum_omega_names,
  multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names,
  page_length_grid_size = page_length_grid_size, weighted_length_grid_size = weighted_length_grid_size,
  finite_eval_grid_size = finite_eval_grid_size, open_eval_grid_size = open_eval_grid_size, batch_size = batch_size, overwrite_existing = TRUE)
write_run_summary(summary_lines, OUTPUT_DIR)
prepare_fresh_alt_output_file(OUTPUT_DIR, backup_existing = backup_existing_outputs)
cv <- empty_alt_table(); start_time <- Sys.time()
finite_T <- T_grid[is.finite(T_grid)]
for (T in finite_T) {
  sim_T <- simulate_finite_alt_critical_values_T(T = T, q_grid = q_grid, gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec,
    weighted_cusum_omega_names = weighted_cusum_omega_names, multiscale_h_sets = multiscale_h_sets,
    multiscale_scale_names = multiscale_scale_names, alpha_levels = alpha_levels, nrep = nrep_finite,
    grid_monitor = grid_monitor, page_length_grid_size = page_length_grid_size,
    weighted_length_grid_size = weighted_length_grid_size, finite_eval_grid_size = finite_eval_grid_size,
    ridge = ridge, range_floor = range_floor, ncores = ncores, seed = seed_finite + as.integer(round(100 * T)),
    exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan, batch_size = batch_size, verbose = verbose)
  cv <- append_alt_rows(cv, sim_T$rows); save_alt_critical_values(cv, OUTPUT_DIR)
}
if (any(!is.finite(T_grid))) {
  sim_inf <- simulate_openend_alt_critical_values(q_grid = q_grid, gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec,
    weighted_cusum_omega_names = weighted_cusum_omega_names, multiscale_h_sets = multiscale_h_sets,
    multiscale_scale_names = multiscale_scale_names, alpha_levels = alpha_levels, nrep = nrep_openend,
    n_train_grid = n_train_grid_open, n_open_grid = n_open_grid_open, open_pair_lag_grid_size = open_pair_lag_grid_size,
    open_eval_grid_size = open_eval_grid_size, ridge = ridge, range_floor = range_floor, ncores = ncores,
    seed = seed_openend, exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan,
    batch_size = batch_size, verbose = verbose)
  cv <- append_alt_rows(cv, sim_inf$rows); save_alt_critical_values(cv, OUTPUT_DIR)
}
paths <- save_alt_critical_values(cv, OUTPUT_DIR)
elapsed <- difftime(Sys.time(), start_time, units = "mins")
message("Done. Freshly simulated alternative-detector critical values written to:")
message("  ", normalize_path2(paths$path, mustWork = FALSE))
message(sprintf("Elapsed time: %.2f minutes", as.numeric(elapsed)))
