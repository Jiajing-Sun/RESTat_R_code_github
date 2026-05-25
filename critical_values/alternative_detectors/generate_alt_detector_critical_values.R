# ==============================================================
# generate_alt_detector_critical_values.R
# ============================================================== 

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
source_project_paths_bootstrap <- function() {
  p <- bootstrap_script_path()
  candidates <- unique(c(
    if (!is.null(p)) file.path(dirname(p), "R", "project_paths.R") else NULL,
    file.path(getwd(), "R", "project_paths.R"),
    file.path(dirname(getwd()), "R", "project_paths.R")
  ))
  for (x in candidates) {
    if (!is.null(x) && file.exists(x)) {
      source(x, local = FALSE)
      return(invisible(x))
    }
  }
  stop("Could not locate R/project_paths.R during bootstrap.")
}
source_project_paths_bootstrap()
ROOT <- resolve_project_root(default_start = dirname(bootstrap_script_path() %||% getwd()))
source(file.path(ROOT, "R", "project_paths.R"), local = FALSE)
source(file.path(ROOT, "R", "dependencies.R"), local = FALSE)
source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "simulation_settings.R"), local = FALSE)
source(file.path(ROOT, "R", "alt_detector_weights.R"), local = FALSE)
source(file.path(ROOT, "R", "method_catalog.R"), local = FALSE)
source(file.path(ROOT, "R", "critical_values_lookup.R"), local = FALSE)
source(file.path(ROOT, "R", "alt_critical_value_generator.R"), local = FALSE)

args <- parse_named_args()
install_missing <- identical(tolower(args$install_missing %||% "false"), "true")
ensure_simulation_packages(install_if_missing = install_missing)

q_grid <- parse_integer_vector_arg(args$q_grid, default = 1:30)
alpha_levels <- parse_numeric_vector_arg(args$alpha_levels, default = c(0.01, 0.05, 0.10))
T_grid <- parse_T_vector_arg(args$T_grid, default = c(1, 2, 5, 10))
include_openend <- parse_bool_arg(args$include_openend, default = FALSE)
if (isTRUE(include_openend) && !any(!is.finite(T_grid))) T_grid <- c(T_grid, Inf)
if (!isTRUE(include_openend)) T_grid <- T_grid[is.finite(T_grid)]
finite_nrep <- parse_scalar_integer_arg(args$nrep_finite %||% args$finite_nrep, default = 10000L)
open_nrep <- parse_scalar_integer_arg(args$nrep_openend %||% args$open_nrep, default = 5000L)
finite_kmax <- parse_scalar_integer_arg(args$finite_kmax, default = 10000L)
open_train_grid_size <- parse_scalar_integer_arg(args$open_train_grid_size, default = 1500L)
open_monitor_grid_size <- parse_scalar_integer_arg(args$open_monitor_grid_size, default = 2000L)
open_s_max <- parse_scalar_numeric_arg(args$open_s_max, default = 20)
page_length_grid_size <- parse_scalar_integer_arg(args$page_length_grid_size, default = 40L)
weighted_length_grid_size <- parse_scalar_integer_arg(args$weighted_length_grid_size, default = 40L)
progress_every <- parse_scalar_integer_arg(args$progress_every, default = 100L)
ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg)
chunk_size <- parse_scalar_integer_arg(args$chunk_size, default = 100L)
overwrite <- parse_bool_arg(args$overwrite, default = TRUE)

output_path <- args$output %||% file.path(ROOT, "outputs", "critical_values_alt_detectors.csv")
ensure_dir(dirname(output_path))
if (file.exists(output_path) && !overwrite) {
  stop("Output file already exists and --overwrite=FALSE: ", output_path)
}
backup_file_if_exists(output_path)

out <- generate_alt_detector_critical_values(
  T_grid = T_grid,
  q_grid = q_grid,
  alpha_levels = alpha_levels,
  finite_nrep = finite_nrep,
  open_nrep = open_nrep,
  finite_kmax = finite_kmax,
  open_train_grid_size = open_train_grid_size,
  open_monitor_grid_size = open_monitor_grid_size,
  open_s_max = open_s_max,
  page_length_grid_size = page_length_grid_size,
  weighted_length_grid_size = weighted_length_grid_size,
  progress_every = progress_every,
  ncores = core_plan$used,
  chunk_size = chunk_size,
  root = ROOT
)
write_csv_atomic(out, output_path, row.names = FALSE)

alt_catalog <- build_alt_method_catalog(
  gamma_vec = c(0, 0.15),
  mosum_h_vec = c(0.10, 0.20),
  weighted_omega_names = c("InvSqrt"),
  multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
  multiscale_scale_names = c("Equal")
)
cv_alt <- normalize_alt_critical_values(read_csv_strict(output_path))
validate_alt_critical_value_coverage(
  cv_alt = cv_alt,
  alt_catalog = alt_catalog,
  T_grid = T_grid,
  q_grid = q_grid,
  alpha_levels = alpha_levels
)

lines <- c(
  strrep("=", 78L),
  "Benchmark-detector critical values regenerated",
  sprintf("Output path: %s", normalize_path2(output_path, mustWork = FALSE)),
  sprintf("Rows written: %d", nrow(out)),
  sprintf("q-grid: %s", format_q_grid(q_grid)),
  sprintf("T-grid: %s", paste(canonical_T_vec(T_grid), collapse = ", ")),
  sprintf("alpha levels: %s", paste(alpha_levels, collapse = ", ")),
  sprintf("Finite/Open nrep: %d / %d", finite_nrep, open_nrep),
  sprintf("Workers requested/used: %d/%d", core_plan$requested, core_plan$used),
  sprintf("Chunk size: %d", chunk_size),
  "Coverage validation: passed",
  strrep("=", 78L)
)
writeLines(lines)
