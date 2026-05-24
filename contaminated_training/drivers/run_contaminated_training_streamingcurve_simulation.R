# ==============================================================
# run_contaminated_training_streamingcurve_simulation.R
# Windows-safe contaminated Phase-I training simulation
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
source(file.path(ROOT, "R", "critical_values_lookup.R"), local = FALSE)
source(file.path(ROOT, "R", "method_catalog.R"), local = FALSE)
source(file.path(ROOT, "R", "simulation_core.R"), local = FALSE)
source(file.path(ROOT, "R", "fpca_pipeline.R"), local = FALSE)
source(file.path(ROOT, "R", "scenarios.R"), local = FALSE)
source(file.path(ROOT, "R", "contamination.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_main.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_alt.R"), local = FALSE)
source(file.path(ROOT, "R", "genData.R"), local = FALSE)

args <- parse_named_args()
install_missing <- identical(tolower(args$install_missing %||% "false"), "true")
ensure_simulation_packages(install_if_missing = install_missing)
output_tier <- sanitize_tag(args$output_tier %||% "final_results")
OUTPUT_ROOT <- file.path(ROOT, "outputs", output_tier)

# ----------------------- user controls -----------------------
# DFT-style representative contamination experiment
scenario <- args$scenario %||% "level_shift"
m_vals <- if (!is.null(args$m)) as.integer(strsplit(args$m, ",")[[1L]]) else c(1000L)
T_grid <- if (!is.null(args$T)) as.numeric(strsplit(args$T, ",")[[1L]]) else c(2)
dgp_types <- if (!is.null(args$dgp)) strsplit(args$dgp, ",")[[1L]] else c("BB", "fIID", "fMA1")
break_fracs <- if (!is.null(args$break_fracs)) as.numeric(strsplit(args$break_fracs, ",")[[1L]]) else c(0.1, 0.8)
train_contam_types <- if (!is.null(args$train_contam)) strsplit(args$train_contam, ",")[[1L]] else c("late_break", "drift")
train_b_grid <- if (!is.null(args$train_b)) as.numeric(strsplit(args$train_b, ",")[[1L]]) else c(0, 0.002, 0.005, 0.010, 0.020, 0.030)
alpha_levels <- if (!is.null(args$alpha_levels)) as.numeric(strsplit(args$alpha_levels, ",")[[1L]]) else c(0.10, 0.05, 0.01)
gamma_vec <- c(0, 0.15)
cvm_weights <- c("U", "Early", "Mid", "Late")
selected_alt_detectors <- if (!is.null(args$selected_alt_detectors)) strsplit(args$selected_alt_detectors, ",")[[1L]] else c("PageCUSUM", "WeightedCUSUM", "MultiscaleMOSUM")
selected_alt_detectors <- trimws(selected_alt_detectors)
selected_alt_detectors <- selected_alt_detectors[nzchar(selected_alt_detectors)]
include_alt_detectors <- length(selected_alt_detectors) > 0L
mosum_h_vec <- c(0.10, 0.20)
weighted_omega_names <- c("InvSqrt")
multiscale_h_sets <- list(H050_100_200 = c(0.05, 0.10, 0.20))
multiscale_scale_names <- c("Equal")

# Representative Phase-II alternative, mirroring the DFT contamination appendix
select_contam_delta <- function(dgp_type, scenario) {
  d <- toupper(as.character(dgp_type))
  sc <- tolower(as.character(scenario))
  if (sc %in% c("level_shift", "smooth_change")) {
    if (d == "BB") return(0.095)
    return(0.017)
  }
  if (sc == "abrupt_local_change") {
    if (d == "BB") return(0.700)
    return(0.110)
  }
  stop("Unknown scenario for contaminated-training driver: ", scenario)
}

nsim <- as.integer(args$nsim %||% 1000L)
basis_k <- 5L
q_cap <- 30L
fve_threshold <- 0.95
fixed_q <- NA_integer_
n_grid <- 301L
nbasis <- 21L
basis_type <- "bspline"
muInfo <- list(type = "sin", a = 1)
varType <- "A"
gaussian <- TRUE
kappa <- 0.7
ridge <- 1e-10
range_floor <- 1e-8
hac_bandwidth <- NULL
page_length_grid_size <- 40L
weighted_length_grid_size <- 40L
finite_eval_grid_size <- 250L
exact_page_scan <- FALSE
exact_weighted_scan <- FALSE
batch_size <- as.integer(args$batch_size %||% 5L)
overwrite_existing <- identical(tolower(args$overwrite %||% "false"), "true")

# Contamination design defaults
train_break_frac_default <- as.numeric(args$train_break_frac %||% 0.8)
train_drift_start_default <- as.numeric(args$train_drift_start_frac %||% 0.0)
contam_shape <- args$contam_shape %||% "phase2_matched"

ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg)
ncores <- core_plan$used

BATCH_ROOT <- file.path(OUTPUT_ROOT, "batches", sanitize_tag(scenario))
RAW_ROOT <- file.path(OUTPUT_ROOT, "raw", sanitize_tag(scenario))
SUMMARY_DIR <- file.path(OUTPUT_ROOT, "summary")
ensure_dir(BATCH_ROOT)
ensure_dir(RAW_ROOT)
ensure_dir(SUMMARY_DIR)

cv_main <- load_main_critical_values(ROOT)
cv_alt <- load_alt_critical_values(ROOT)
main_catalog <- build_main_method_catalog(gamma_vec = gamma_vec, cvm_weights = cvm_weights)
alt_catalog <- if (include_alt_detectors) {
  build_alt_method_catalog(gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec,
                           weighted_omega_names = weighted_omega_names,
                           multiscale_h_sets = multiscale_h_sets,
                           multiscale_scale_names = multiscale_scale_names)
} else {
  data.frame()
}
if (nrow(alt_catalog) > 0L) {
  alt_catalog <- subset(alt_catalog, detector %in% selected_alt_detectors)
}

summary_lines <- c(
  strrep("=", 78L),
  "Streaming-curve contaminated-training simulation (Windows-safe build)",
  sprintf("Project root: %s", normalize_path2(ROOT, mustWork = FALSE)),
  sprintf("Scenario: %s", scenario),
  sprintf("Output tier: %s", output_tier),
  sprintf("DGP types: %s", paste(dgp_types, collapse = ", ")),
  sprintf("m-grid: %s", paste(m_vals, collapse = ", ")),
  sprintf("T-grid: %s", paste(T_grid, collapse = ", ")),
  sprintf("break fractions: %s", paste(break_fracs, collapse = ", ")),
  sprintf("training contamination types: %s", paste(train_contam_types, collapse = ", ")),
  sprintf("train_b grid: %s", paste(train_b_grid, collapse = ", ")),
  sprintf("alpha levels: %s", paste(alpha_levels, collapse = ", ")),
  sprintf("retained alternative detectors: %s", if (length(selected_alt_detectors) > 0L) paste(selected_alt_detectors, collapse = ", ") else "<none>"),
  sprintf("Replications per parameter combination: %d", nsim),
  sprintf("Workers requested/used: %d/%d", core_plan$requested, core_plan$used),
  sprintf("Windows-safe connection budget: %d", core_plan$connection_budget),
  sprintf("Windows-safe hard cap: %d", core_plan$os_cap),
  strrep("=", 78L)
)
write_run_summary(summary_lines, SUMMARY_DIR, filename = paste0("run_summary_train_contam_", sanitize_tag(scenario), ".txt"))

run_chunk <- function(seed_chunk, sim_par, cv_main, cv_alt, main_catalog, alt_catalog) {
  out <- lapply(seed_chunk, function(sd) {
    set.seed(sd)
    rep_id <- as.integer(sd)
    simulate_one_replicate(rep_id = rep_id, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                           main_catalog = main_catalog, alt_catalog = alt_catalog)
  })
  do.call(rbind, out)
}

param_grid <- expand.grid(m = m_vals,
                          T_value = T_grid,
                          dgp_type = dgp_types,
                          break_frac = break_fracs,
                          train_contam = train_contam_types,
                          train_b = train_b_grid,
                          stringsAsFactors = FALSE)
param_grid <- param_grid[order(param_grid$dgp_type, param_grid$m, param_grid$T_value,
                               param_grid$break_frac, param_grid$train_contam, param_grid$train_b), , drop = FALSE]
row.names(param_grid) <- NULL

for (i in seq_len(nrow(param_grid))) {
  combo <- param_grid[i, , drop = FALSE]
  delta_use <- select_contam_delta(combo$dgp_type, scenario)
  s_star_use <- as.integer(round(as.numeric(combo$break_frac) * as.numeric(combo$m) * as.numeric(combo$T_value)))
  file_tag <- sanitize_tag(sprintf(
    "traincontam_%s_%s_m%s_T%s_bf%s_ct%s_tb%s_nsim%s",
    scenario, combo$dgp_type, combo$m, canonical_T_scalar(combo$T_value),
    format(combo$break_frac, trim = TRUE), combo$train_contam, format(combo$train_b, trim = TRUE), nsim
  ))
  batch_dir <- file.path(BATCH_ROOT, file_tag)
  final_rds <- file.path(RAW_ROOT, paste0(file_tag, ".rds"))
  final_csv <- file.path(RAW_ROOT, paste0(file_tag, ".csv"))
  progress_log <- file.path(batch_dir, "progress.log")
  ensure_dir(batch_dir)

  sim_par <- list(
    scenario = scenario,
    dgp_type = as.character(combo$dgp_type),
    m = as.integer(combo$m),
    T_value = as.numeric(combo$T_value),
    delta = delta_use,
    s_star = s_star_use,
    basis_k = basis_k,
    alpha_levels = alpha_levels,
    q_cap = q_cap,
    fve_threshold = fve_threshold,
    fixed_q = fixed_q,
    n_grid = n_grid,
    nbasis = nbasis,
    basis_type = basis_type,
    muInfo = muInfo,
    varType = varType,
    gaussian = gaussian,
    kappa = kappa,
    ridge = ridge,
    range_floor = range_floor,
    hac_bandwidth = hac_bandwidth,
    page_length_grid_size = page_length_grid_size,
    weighted_length_grid_size = weighted_length_grid_size,
    finite_eval_grid_size = finite_eval_grid_size,
    exact_page_scan = exact_page_scan,
    exact_weighted_scan = exact_weighted_scan,
    multiscale_h_sets = multiscale_h_sets,
    train_contam = as.character(combo$train_contam),
    train_b = as.numeric(combo$train_b),
    train_break_frac = if (as.character(combo$train_contam) == "late_break") train_break_frac_default else NA_real_,
    train_drift_start_frac = if (as.character(combo$train_contam) == "drift") train_drift_start_default else NA_real_,
    contam_shape = contam_shape,
    t_train_star = NA_integer_
  )

  message(sprintf("Contaminated-training simulation: DGP=%s | m=%d | T=%s | break_frac=%.2f | contam=%s | b_train=%.3f | delta=%.3f | nsim=%d | workers=%d",
                  sim_par$dgp_type, sim_par$m, canonical_T_scalar(sim_par$T_value), as.numeric(combo$break_frac),
                  sim_par$train_contam, sim_par$train_b, sim_par$delta, nsim, ncores))

  seeds <- make_seed_stream(nsim, seed = 20260403L + i * 10000L)
  seed_chunks <- split(seeds, ceiling(seq_along(seeds) / batch_size))
  batch_files <- file.path(batch_dir, sprintf("batch_%04d.rds", seq_along(seed_chunks)))

  cl <- NULL
  if (ncores > 1L) {
    cl <- make_psock_cluster(min(ncores, batch_size))
    on.exit(stop_psock_cluster(cl), add = TRUE)
    ROOT2 <- ROOT
    parallel::clusterExport(cl, varlist = c("ROOT2", "sim_par", "cv_main", "cv_alt", "main_catalog", "alt_catalog", "run_chunk"), envir = environment())
    parallel::clusterEvalQ(cl, {
      source(file.path(ROOT2, "R", "worker_bundle.R"), local = FALSE)
      NULL
    })
  }

  for (b in seq_along(seed_chunks)) {
    if (file_nonempty(batch_files[b]) && !overwrite_existing) {
      append_progress_line(progress_log, sprintf("skip existing batch %04d", b))
      next
    }
    res_batch <- if (is.null(cl)) {
      run_chunk(seed_chunks[[b]], sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                main_catalog = main_catalog, alt_catalog = alt_catalog)
    } else {
      parts <- split(seed_chunks[[b]], ceiling(seq_along(seed_chunks[[b]]) / max(1L, ceiling(length(seed_chunks[[b]]) / ncores))))
      res_list <- parallel::parLapplyLB(cl, parts, function(sd_chunk) {
        run_chunk(sd_chunk, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                  main_catalog = main_catalog, alt_catalog = alt_catalog)
      })
      combine_rbind(res_list)
    }
    atomic_save_rds(res_batch, batch_files[b])
    append_progress_line(progress_log, sprintf("completed batch %04d of %04d (%d rows)", b, length(seed_chunks), nrow(res_batch)))
  }

  stop_psock_cluster(cl)
  on.exit(NULL, add = FALSE)

  batch_objs <- lapply(batch_files[file_nonempty(batch_files)], safe_read_rds)
  final_df <- combine_rbind(batch_objs)
  atomic_save_rds(final_df, final_rds)
  atomic_write_csv(final_df, final_csv)
}

message("Contaminated-training simulation complete. Batch files are in: ", normalize_path2(BATCH_ROOT, mustWork = FALSE))
message("Combined raw files are in: ", normalize_path2(RAW_ROOT, mustWork = FALSE))
