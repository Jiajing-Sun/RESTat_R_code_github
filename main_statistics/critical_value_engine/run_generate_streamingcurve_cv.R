# ==============================================================
# run_generate_streamingcurve_cv.R
# ============================================================== 
# Fresh simulation driver for the streaming-curve paper.
# No warm starts, no top-up logic, no use of older critical-value tables.
# Includes SSMS / RSMS / HAC KS and weighted-CvM critical values.

bootstrap_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) {
    return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  }

  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p1)) return(normalizePath(p1, winslash = "/", mustWork = FALSE))
    p2 <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }

  NULL
}

source_bootstrap_file <- function(rel_path) {
  candidates <- unique(c(
    file.path(dirname(bootstrap_script_path() %||% getwd()), rel_path),
    file.path(getwd(), rel_path)
  ))

  for (p in candidates) {
    if (file.exists(p)) {
      source(p, local = FALSE)
      return(invisible(p))
    }
  }

  stop("Could not source bootstrap helper: ", rel_path)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
source_bootstrap_file(file.path("R", "project_paths.R"))

BOOTSTRAP_PATH <- bootstrap_script_path()
BOOTSTRAP_DIR <- if (!is.null(BOOTSTRAP_PATH)) dirname(BOOTSTRAP_PATH) else getwd()
ROOT <- resolve_project_root(default_start = BOOTSTRAP_DIR)

source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "weights.R"), local = FALSE)
source(file.path(ROOT, "R", "finite_critical_values.R"), local = FALSE)
source(file.path(ROOT, "R", "openend_critical_values.R"), local = FALSE)
source(file.path(ROOT, "R", "critical_values_io.R"), local = FALSE)

# --------------------------------------------------------------
# User controls
# --------------------------------------------------------------
q_grid <- 1:30
T_grid <- c(1, 2, 5, 10, Inf)
gamma_vec <- c(0, 0.15)
weight_names <- c("U", "Early", "Mid", "Late")
alpha_levels <- c(0.10, 0.05, 0.01)

# Monte Carlo controls
nrep_finite <- 10000L
nrep_openend <- 5000L

grid_monitor <- 10000L      # finite T: k_max = m * T = grid_monitor
n_train_grid_open <- 1500L  # open-end training bridge grid
n_open_grid_open <- 2000L   # open-end x-grid on [0,1]

ridge <- 1e-10
range_floor <- 1e-8

args <- parse_named_args()
ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg)
ncores <- core_plan$used

seed_finite <- 123456789L
seed_openend <- 13579L
verbose <- TRUE
backup_existing_outputs <- TRUE

# Paths
OUTPUT_DIR <- file.path(ROOT, "outputs")
BASE_OUT <- file.path(OUTPUT_DIR, "critical_values_all.csv")
WEIGHTS_OUT <- file.path(OUTPUT_DIR, "critical_values_all_weights.csv")
ensure_dir(OUTPUT_DIR)

# --------------------------------------------------------------
# Fresh-run reporting and output preparation
# --------------------------------------------------------------
summary_lines <- report_simulation_design(
  root = ROOT,
  output_dir = OUTPUT_DIR,
  q_grid = q_grid,
  T_grid = T_grid,
  gamma_vec = gamma_vec,
  weight_names = weight_names,
  alpha_levels = alpha_levels,
  nrep_finite = nrep_finite,
  nrep_openend = nrep_openend,
  grid_monitor = grid_monitor,
  n_train_grid_open = n_train_grid_open,
  n_open_grid_open = n_open_grid_open,
  core_plan = core_plan,
  overwrite_existing = TRUE
)
write_run_summary(summary_lines, OUTPUT_DIR)
prepare_fresh_output_files(OUTPUT_DIR, backup_existing = backup_existing_outputs)

# --------------------------------------------------------------
# Simulate everything from scratch
# --------------------------------------------------------------
cv <- list(base = empty_base_table(), weights = empty_weight_table())
start_time <- Sys.time()

finite_T <- T_grid[is.finite(T_grid)]
for (T in finite_T) {
  sim_T <- simulate_finite_critical_values_T(
    T = T,
    q_grid = q_grid,
    gamma_vec = gamma_vec,
    weight_names = weight_names,
    alpha_levels = alpha_levels,
    nrep = nrep_finite,
    grid_monitor = grid_monitor,
    ridge = ridge,
    range_floor = range_floor,
    ncores = ncores,
    seed = seed_finite + as.integer(round(100 * T)),
    verbose = verbose
  )

  cv <- append_cv_rows(cv, sim_T$base_rows, sim_T$weight_rows)
  save_critical_values(cv, OUTPUT_DIR)
}

if (any(!is.finite(T_grid))) {
  sim_inf <- simulate_openend_critical_values(
    q_grid = q_grid,
    gamma_vec = gamma_vec,
    weight_names = weight_names,
    alpha_levels = alpha_levels,
    nrep = nrep_openend,
    n_train_grid = n_train_grid_open,
    n_open_grid = n_open_grid_open,
    ridge = ridge,
    range_floor = range_floor,
    ncores = ncores,
    seed = seed_openend,
    verbose = verbose
  )

  cv <- append_cv_rows(cv, sim_inf$base_rows, sim_inf$weight_rows)
  save_critical_values(cv, OUTPUT_DIR)
}

# --------------------------------------------------------------
# Final summary
# --------------------------------------------------------------
paths <- save_critical_values(cv, OUTPUT_DIR)
elapsed <- difftime(Sys.time(), start_time, units = "mins")

message("Done. Freshly simulated critical-value tables written to:")
message("  ", normalize_path2(paths$base_path, mustWork = FALSE))
message("  ", normalize_path2(paths$weights_path, mustWork = FALSE))
message(sprintf("Elapsed time: %.2f minutes", as.numeric(elapsed)))
