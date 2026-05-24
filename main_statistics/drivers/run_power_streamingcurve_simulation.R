# ==============================================================
# run_power_streamingcurve_simulation.R
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
source(file.path(ROOT, "R", "critical_values_lookup.R"), local = FALSE)
source(file.path(ROOT, "R", "method_catalog.R"), local = FALSE)
source(file.path(ROOT, "R", "preflight_checks.R"), local = FALSE)
source(file.path(ROOT, "R", "simulation_core.R"), local = FALSE)

args <- parse_named_args()
install_missing <- identical(tolower(args$install_missing %||% "false"), "true")
ensure_simulation_packages(install_if_missing = install_missing)
source(file.path(ROOT, "R", "fpca_pipeline.R"), local = FALSE)
source(file.path(ROOT, "R", "scenarios.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_main.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_alt.R"), local = FALSE)
source(file.path(ROOT, "R", "genData.R"), local = FALSE)

# ----------------------- user controls -----------------------
settings <- resolve_power_run_settings(args)

ncores_arg <- if (!is.null(args$ncores)) as.integer(args$ncores) else NULL
core_plan <- resolve_safe_ncores(requested = ncores_arg)
ncores <- core_plan$used

OUTPUT_DIR <- file.path(ROOT, "outputs", "power_raw", sanitize_tag(settings$scenario))
ensure_dir(OUTPUT_DIR)

main_catalog <- build_main_method_catalog(gamma_vec = settings$gamma_vec, cvm_weights = settings$cvm_weights)
alt_catalog <- if (settings$include_alt_detectors) {
  build_alt_method_catalog(gamma_vec = settings$gamma_vec, mosum_h_vec = settings$mosum_h_vec,
                           weighted_omega_names = settings$weighted_omega_names,
                           multiscale_h_sets = settings$multiscale_h_sets,
                           multiscale_scale_names = settings$multiscale_scale_names)
} else {
  data.frame()
}
preflight <- run_simulation_preflight(
  root = ROOT,
  settings = settings,
  main_catalog = main_catalog,
  alt_catalog = alt_catalog,
  ncores = ncores,
  check_worker_bootstrap = ncores > 1L
)
cv_main <- preflight$cv_main
cv_alt <- preflight$cv_alt
if (is.null(cv_alt)) alt_catalog <- data.frame()
if (!isTRUE(preflight$worker_bootstrap_ok) && ncores > 1L) {
  message("PSOCK bootstrap check failed; falling back to serial execution for this run.")
  ncores <- 1L
}

summary_lines <- c(
  strrep("=", 78L),
  sprintf("Streaming-curve power simulation [%s]", settings$scenario),
  sprintf("Project root: %s", normalize_path2(ROOT, mustWork = FALSE)),
  sprintf("Output directory: %s", normalize_path2(OUTPUT_DIR, mustWork = FALSE)),
  sprintf("m-grid: %s", paste(settings$m_vals, collapse = ", ")),
  sprintf("T-grid: %s", paste(settings$T_grid, collapse = ", ")),
  sprintf("DGP types: %s", paste(settings$dgp_types, collapse = ", ")),
  sprintf("s* grid: %s", paste(settings$s_star_vals, collapse = ", ")),
  sprintf("delta grid by DGP: %s", format_delta_map_summary(settings$delta_map)),
  sprintf("Include benchmark detectors: %s", nrow(alt_catalog) > 0L),
  sprintf("alpha levels: %s", paste(settings$alpha_levels, collapse = ", ")),
  sprintf("Replications per parameter combination: %d", settings$nsim),
  sprintf("Workers requested/used: %d/%d", core_plan$requested, ncores),
  "Preflight checks: passed",
  strrep("=", 78L)
)
if (!is.null(preflight$alt_warning)) {
  summary_lines <- append(summary_lines, sprintf("Benchmark detector status: %s", preflight$alt_warning), after = 10L)
}
if (!isTRUE(preflight$worker_bootstrap_ok) && core_plan$used > 1L) {
  summary_lines <- append(summary_lines, "Parallel status: PSOCK unavailable during preflight, using serial fallback.", after = 11L)
}
write_run_summary(summary_lines, OUTPUT_DIR, filename = "run_summary_power.txt")

run_chunk <- function(seed_chunk, sim_par, cv_main, cv_alt, main_catalog, alt_catalog) {
  out <- lapply(seed_chunk, function(sd) {
    set.seed(sd)
    rep_id <- as.integer(sd)
    simulate_one_replicate(rep_id = rep_id, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                           main_catalog = main_catalog, alt_catalog = alt_catalog)
  })
  do.call(rbind, out)
}

param_grid <- build_power_param_grid(settings)
for (i in seq_len(nrow(param_grid))) {
  combo <- param_grid[i, , drop = FALSE]
  file_tag <- sanitize_tag(sprintf("%s_%s_m%s_T%s_delta%s_sstar%s_nsim%s",
                                   settings$scenario, combo$dgp_type, combo$m, canonical_T_scalar(combo$T_value),
                                   combo$delta, combo$s_star, settings$nsim))
  out_file <- file.path(OUTPUT_DIR, paste0(file_tag, ".csv"))
  if (file.exists(out_file) && !settings$overwrite_existing) {
    message("Skipping existing file: ", basename(out_file))
    next
  }

  sim_par <- list(
    scenario = settings$scenario,
    dgp_type = as.character(combo$dgp_type),
    m = as.integer(combo$m),
    T_value = as.numeric(combo$T_value),
    delta = as.numeric(combo$delta),
    s_star = as.integer(combo$s_star),
    basis_k = settings$basis_k,
    alpha_levels = settings$alpha_levels,
    q_cap = settings$q_cap,
    fve_threshold = settings$fve_threshold,
    fixed_q = settings$fixed_q,
    n_grid = settings$n_grid,
    ridge = settings$ridge,
    range_floor = settings$range_floor,
    hac_bandwidth = settings$hac_bandwidth,
    page_length_grid_size = settings$page_length_grid_size,
    weighted_length_grid_size = settings$weighted_length_grid_size,
    finite_eval_grid_size = settings$finite_eval_grid_size,
    exact_page_scan = settings$exact_page_scan,
    exact_weighted_scan = settings$exact_weighted_scan,
    multiscale_h_sets = settings$multiscale_h_sets
  )

  message(sprintf("Running power simulation: scenario=%s | DGP=%s | m=%d | T=%s | delta=%s | s*=%d | nsim=%d | workers=%d",
                  sim_par$scenario, sim_par$dgp_type, sim_par$m, canonical_T_scalar(sim_par$T_value),
                  format(sim_par$delta, trim = TRUE), sim_par$s_star, settings$nsim, ncores))

  seeds <- make_seed_stream(settings$nsim, seed = 20260331L + i * 10000L)
  seed_chunks <- split(seeds, ceiling(seq_along(seeds) / settings$batch_size))

  if (ncores <= 1L || length(seed_chunks) == 1L) {
    res_list <- lapply(seed_chunks, run_chunk, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                       main_catalog = main_catalog, alt_catalog = alt_catalog)
  } else {
    cluster_error <- NULL
    cl <- tryCatch(
      make_psock_cluster(min(ncores, length(seed_chunks))),
      error = function(e) {
        cluster_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(cl)) {
      message("PSOCK cluster unavailable for this parameter combination; falling back to serial execution. Details: ", cluster_error)
      res_list <- lapply(seed_chunks, run_chunk, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                         main_catalog = main_catalog, alt_catalog = alt_catalog)
    } else {
      on.exit(stop_psock_cluster(cl), add = TRUE)
      ROOT2 <- ROOT
      parallel::clusterExport(cl, varlist = c("ROOT2", "sim_par", "cv_main", "cv_alt", "main_catalog", "alt_catalog", "run_chunk"), envir = environment())
      parallel::clusterEvalQ(cl, {
        source(file.path(ROOT2, "R", "worker_bundle.R"), local = FALSE)
        NULL
      })
      res_list <- parallel::parLapplyLB(cl, seed_chunks, function(sd_chunk) {
        run_chunk(sd_chunk, sim_par = sim_par, cv_main = cv_main, cv_alt = cv_alt,
                  main_catalog = main_catalog, alt_catalog = alt_catalog)
      })
      stop_psock_cluster(cl)
      on.exit(NULL, add = FALSE)
    }
  }

  res <- do.call(rbind, res_list)
  write_csv_atomic(res, out_file, row.names = FALSE)
}

message("Power simulation complete. Raw files are in: ", normalize_path2(OUTPUT_DIR, mustWork = FALSE))
