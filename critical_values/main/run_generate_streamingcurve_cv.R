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

parse_csv_tokens <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(character(0L))
  tokens <- unlist(strsplit(as.character(x[1L]), ",", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens[nzchar(tokens)]
}

parse_bool_arg <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(isTRUE(default))
  key <- toupper(trimws(as.character(x[1L])))
  if (key %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (key %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop("Invalid logical value: ", x[1L])
}

parse_character_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  unique(tokens)
}

parse_numeric_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- suppressWarnings(as.numeric(tokens))
  if (any(is.na(vals))) stop("Invalid numeric vector argument: ", as.character(x[1L]))
  unique(vals)
}

parse_integer_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- unlist(lapply(tokens, function(tok) {
    if (grepl("^-?[0-9]+:-?[0-9]+$", tok)) {
      endpoints <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1L]])
      return(seq(endpoints[1L], endpoints[2L]))
    }
    val <- suppressWarnings(as.numeric(tok))
    if (is.na(val) || abs(val - round(val)) > 1e-10) stop("Invalid integer vector argument: ", tok)
    as.integer(round(val))
  }), use.names = FALSE)
  unique(as.integer(vals))
}

parse_T_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- vapply(tokens, function(tok) {
    key <- tolower(trimws(tok))
    if (key %in% c("inf", "infinity")) return(Inf)
    val <- suppressWarnings(as.numeric(tok))
    if (is.na(val)) stop("Invalid T-grid entry: ", tok)
    val
  }, numeric(1L))
  unique(vals)
}

parse_scalar_integer_arg <- function(x, default) {
  vals <- parse_integer_vector_arg(x, default = default)
  as.integer(vals[1L])
}

parse_scalar_numeric_arg <- function(x, default) {
  vals <- parse_numeric_vector_arg(x, default = default)
  as.numeric(vals[1L])
}

# --------------------------------------------------------------
# User controls
# --------------------------------------------------------------
args <- parse_named_args()

q_grid <- parse_integer_vector_arg(args$q_grid, default = 1:30)
if (!is.null(args$q_max)) q_grid <- seq_len(parse_scalar_integer_arg(args$q_max, default = max(q_grid)))
T_grid <- parse_T_vector_arg(args$T_grid, default = c(1, 2, 5, 10))
include_openend <- parse_bool_arg(args$include_openend, default = FALSE)
if (isTRUE(include_openend) && !any(!is.finite(T_grid))) T_grid <- c(T_grid, Inf)
if (!isTRUE(include_openend)) T_grid <- T_grid[is.finite(T_grid)]
gamma_vec <- parse_numeric_vector_arg(args$gamma_vec, default = c(0, 0.15))
weight_names <- parse_character_vector_arg(args$weight_names, default = c("U", "Early", "Mid", "Late"))
alpha_levels <- parse_numeric_vector_arg(args$alpha_levels, default = c(0.10, 0.05, 0.01))

# Monte Carlo controls
nrep_finite <- parse_scalar_integer_arg(args$nrep_finite %||% args$finite_nrep, default = 10000L)
nrep_openend <- parse_scalar_integer_arg(args$nrep_openend %||% args$open_nrep, default = 5000L)

grid_monitor <- parse_scalar_integer_arg(args$grid_monitor %||% args$finite_kmax, default = 10000L)
n_train_grid_open <- parse_scalar_integer_arg(args$n_train_grid_open %||% args$open_train_grid_size, default = 1500L)
n_open_grid_open <- parse_scalar_integer_arg(args$n_open_grid_open %||% args$open_monitor_grid_size, default = 2000L)

ridge <- parse_scalar_numeric_arg(args$ridge, default = 1e-10)
range_floor <- parse_scalar_numeric_arg(args$range_floor, default = 1e-8)

ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg)
ncores <- core_plan$used

seed_finite <- parse_scalar_integer_arg(args$seed_finite, default = 123456789L)
seed_openend <- parse_scalar_integer_arg(args$seed_openend, default = 13579L)
verbose <- parse_bool_arg(args$verbose, default = TRUE)
backup_existing_outputs <- parse_bool_arg(args$backup_existing_outputs, default = TRUE)

# Paths
OUTPUT_DIR <- args$output_dir %||% file.path(ROOT, "outputs")
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
