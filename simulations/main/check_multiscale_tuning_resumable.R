# ==============================================================
# check_multiscale_tuning_resumable.R
# Resumable tuning diagnostic for multiscale MOSUM sensitivity
# ==============================================================

bootstrap_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  NULL
}

source_project_paths_bootstrap <- function() {
  p <- bootstrap_script_path()
  candidates <- unique(c(
    if (!is.null(p)) file.path(dirname(p), "R", "project_paths.R") else NULL,
    if (!is.null(p)) file.path(dirname(p), "..", "R", "project_paths.R") else NULL,
    file.path(getwd(), "R", "project_paths.R"),
    file.path(dirname(getwd()), "R", "project_paths.R"),
    file.path(getwd(), "..", "R", "project_paths.R")
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
SCRIPT_PATH <- bootstrap_script_path() %||% file.path(getwd(), "drivers", "check_multiscale_tuning_targeted.R")
ROOT <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_paths.R"), local = FALSE)
source(file.path(ROOT, "R", "dependencies.R"), local = FALSE)
source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "simulation_settings.R"), local = FALSE)
source(file.path(ROOT, "R", "critical_values_lookup.R"), local = FALSE)
source(file.path(ROOT, "R", "method_catalog.R"), local = FALSE)
source(file.path(ROOT, "R", "simulation_core.R"), local = FALSE)
source(file.path(ROOT, "R", "alt_critical_value_generator.R"), local = FALSE)

ensure_simulation_packages(install_if_missing = FALSE)
source(file.path(ROOT, "R", "fpca_pipeline.R"), local = FALSE)
source(file.path(ROOT, "R", "scenarios.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_main.R"), local = FALSE)
source(file.path(ROOT, "R", "detectors_alt.R"), local = FALSE)
source(file.path(ROOT, "R", "genData.R"), local = FALSE)

args <- parse_named_args()

nsim <- as.integer(args$nsim %||% 1000L)
cv_nrep <- as.integer(args$cv_nrep %||% 5000L)
ncores_req <- as.integer(args$ncores %||% 6L)
batch_size <- as.integer(args$batch_size %||% 20L)
q_cap <- as.integer(args$q_cap %||% 30L)
scenario_vec <- parse_character_vector_arg(args$scenarios, default = c("level_shift", "smooth_change", "abrupt_local_change"))
s_star_vec <- parse_integer_vector_arg(args$s_star_vals, default = c(50L, 200L))
delta_vec_arg <- parse_numeric_vector_arg(args$delta_vals, default = NULL)
standardizer_vec <- parse_character_vector_arg(args$standardizers, default = c("HAC", "SSMS", "RSMS"))
scale_weight_vec <- parse_character_vector_arg(args$scale_weights, default = c("Equal", "InvSqrtH", "SqrtH", "H"))
hset_vec <- parse_character_vector_arg(args$hsets, default = NULL)
open_cv_nrep <- as.integer(args$open_cv_nrep %||% 2000L)
out_tag <- sanitize_tag(args$out_tag %||% sprintf("explore_nsim%s_cv%s", nsim, cv_nrep))
force_cv <- tolower(args$force_cv %||% "false") %in% c("true", "t", "1", "yes")
force_null <- tolower(args$force_null %||% "false") %in% c("true", "t", "1", "yes")

core_plan <- resolve_safe_ncores(requested = ncores_req)
ncores <- core_plan$used

out_dir <- file.path(ROOT, "outputs", "tuning_diagnostics", paste0("multiscale_", out_tag))
combo_dir <- file.path(out_dir, "combo_summaries")
ensure_dir(out_dir)
ensure_dir(combo_dir)

# Override the main-bundle multiscale weighting so the diagnostic can test
# InvSqrtH without editing the production simulation modules.
normalize_scale_weight_name <- function(x) {
  w <- toupper(trimws(as.character(x)[1L]))
  if (w %in% c("EQUAL", "UNIFORM", "ONE", "CONST", "CONSTANT", "1")) return("Equal")
  if (w %in% c("INVSQRTH", "INV_SQRT_H", "1/SQRT_H")) return("InvSqrtH")
  if (w %in% c("INVH", "INV_H", "1/H")) return("InvH")
  if (w %in% c("SQRTH", "SQRT_H", "SQRT(H)")) return("SqrtH")
  if (w %in% c("H", "LINEARH", "LINEAR_H")) return("H")
  if (w %in% c("H2", "H^2", "HSQ", "H_SQ")) return("H2")
  stop("Unknown multiscale weight. Supported: Equal, InvSqrtH, InvH, SqrtH, H, H2")
}

make_multiscale_weight <- function(h, name = "Equal") {
  nm <- normalize_scale_weight_name(name)
  h <- pmax(as.numeric(h), 1e-8)
  if (nm == "Equal") return(rep(1, length(h)))
  if (nm == "InvSqrtH") return(h^(-1/2))
  if (nm == "InvH") return(h^(-1))
  if (nm == "SqrtH") return(h^(1/2))
  if (nm == "H") return(h)
  if (nm == "H2") return(h^2)
  stop("Unknown multiscale weight.")
}

multiscale_h_sets <- list(
  H025_050_075_100 = c(0.025, 0.05, 0.075, 0.10),
  H025_050_100_200 = c(0.025, 0.05, 0.10, 0.20),
  H050_100_200 = c(0.05, 0.10, 0.20),
  H050_075_100_150 = c(0.05, 0.075, 0.10, 0.15),
  H050_100_200_300 = c(0.05, 0.10, 0.20, 0.30),
  H100_150_200_250_300 = c(0.10, 0.15, 0.20, 0.25, 0.30),
  H100_200_300_400 = c(0.10, 0.20, 0.30, 0.40),
  H150_200_250_300 = c(0.15, 0.20, 0.25, 0.30),
  H200_250_300_350_400 = c(0.20, 0.25, 0.30, 0.35, 0.40)
)
if (!is.null(hset_vec)) {
  unknown_hsets <- setdiff(hset_vec, names(multiscale_h_sets))
  if (length(unknown_hsets) > 0L) {
    stop("Unknown hsets: ", paste(unknown_hsets, collapse = ", "),
         ". Available hsets are: ", paste(names(multiscale_h_sets), collapse = ", "))
  }
  multiscale_h_sets <- multiscale_h_sets[hset_vec]
}
multiscale_scale_names <- vapply(scale_weight_vec, normalize_scale_weight_name, character(1L))

alt_catalog <- build_alt_method_catalog(
  gamma_vec = c(0),
  standardizers = standardizer_vec,
  mosum_h_vec = numeric(0),
  weighted_omega_names = character(0),
  multiscale_h_sets = multiscale_h_sets,
  multiscale_scale_names = multiscale_scale_names
)
alt_catalog <- subset(alt_catalog, detector == "MultiscaleMOSUM")

q_grid <- seq_len(q_cap)
alpha_levels <- 0.05

cv_file <- file.path(out_dir, "cv_alt.csv")
if (file.exists(cv_file) && !force_cv) {
  message("Reusing targeted multiscale MOSUM critical values: ", cv_file)
  cv_alt <- read.csv(cv_file, stringsAsFactors = FALSE, check.names = FALSE)
  cv_alt <- normalize_alt_critical_values(cv_alt)
} else {
  message("Generating targeted multiscale MOSUM critical values ...")
  designs <- build_alt_cv_designs(
    T_grid = c(2),
    q_grid = q_grid,
    finite_nrep = cv_nrep,
    open_nrep = open_cv_nrep
  )
  cv_list <- lapply(designs, function(design) {
    sim <- simulate_alt_stat_matrix(
      design = design,
      alt_catalog = alt_catalog,
      q_grid = q_grid,
      ridge = 1e-10,
      range_floor = 1e-8,
      page_length_grid_size = 40L,
      weighted_length_grid_size = 40L,
      exact_page_scan = FALSE,
      exact_weighted_scan = FALSE,
      multiscale_h_sets = multiscale_h_sets,
      progress_every = 500L,
      ncores = 1L,
      chunk_size = 100L,
      root = ROOT
    )
    summarize_alt_stat_matrix(
      stat_matrix = sim$stats,
      column_spec = sim$column_spec,
      alpha_levels = alpha_levels,
      T_label = design$T_label,
      source_file = design$source_file
    )
  })
  cv_alt <- normalize_alt_critical_values(do.call(rbind, cv_list))
  write_csv_atomic(cv_alt, cv_file, row.names = FALSE)
}

run_chunk <- function(seed_chunk, sim_par, cv_alt, alt_catalog) {
  out <- lapply(seed_chunk, function(sd) {
    set.seed(sd)
    simulate_one_replicate(
      rep_id = as.integer(sd),
      sim_par = sim_par,
      cv_main = NULL,
      cv_alt = cv_alt,
      main_catalog = data.frame(),
      alt_catalog = alt_catalog
    )
  })
  do.call(rbind, out)
}

run_combo <- function(sim_par, nsim, ncores, batch_size, cv_alt, alt_catalog, root) {
  seeds <- make_seed_stream(nsim, seed = 20260426L + as.integer(sim_par$m) + 100L * as.integer(sim_par$s_star %||% 0L) + 10L * as.integer(sim_par$delta * 1000))
  seed_chunks <- split(seeds, ceiling(seq_along(seeds) / batch_size))

  if (ncores <= 1L || length(seed_chunks) == 1L) {
    pieces <- lapply(seed_chunks, run_chunk, sim_par = sim_par, cv_alt = cv_alt, alt_catalog = alt_catalog)
    return(do.call(rbind, pieces))
  }
  if (!identical(Sys.info()[["sysname"]], "Windows")) {
    pieces <- parallel::mclapply(
      seed_chunks,
      run_chunk,
      sim_par = sim_par,
      cv_alt = cv_alt,
      alt_catalog = alt_catalog,
      mc.cores = min(ncores, length(seed_chunks))
    )
    return(do.call(rbind, pieces))
  }

  pieces <- lapply(seed_chunks, run_chunk, sim_par = sim_par, cv_alt = cv_alt, alt_catalog = alt_catalog)
  do.call(rbind, pieces)
}

summarize_null <- function(df) {
  ok <- subset(df, (is.na(error_flag) | !error_flag) & !is.na(statistic) & !is.na(reject))
  stats_by_method <- split(ok$statistic, ok$method_id)
  cv_emp <- data.frame(
    method_id = names(stats_by_method),
    empirical_cv = vapply(stats_by_method, function(x) as.numeric(stats::quantile(x, probs = 0.95, names = FALSE, type = 7, na.rm = TRUE)), numeric(1L)),
    stringsAsFactors = FALSE
  )
  size_df <- aggregate(as.numeric(reject) ~ method_id + method_label + standardizer + hset_name + scale_weight_name,
                       data = ok, FUN = mean)
  names(size_df)[names(size_df) == "as.numeric(reject)"] <- "size"
  merge(size_df, cv_emp, by = "method_id", all.x = TRUE, sort = FALSE)
}

summarize_power_combo <- function(df, cv_emp) {
  ok <- subset(df, (is.na(error_flag) | !error_flag) & !is.na(statistic) & !is.na(reject))
  ok <- merge(ok, cv_emp[, c("method_id", "empirical_cv")], by = "method_id", all.x = TRUE, sort = FALSE)
  ok$size_adjusted_reject <- as.numeric(ok$statistic > ok$empirical_cv)
  ok$delay <- ok$first_rejection - ok$s_star
  group_vars <- c("scenario", "dgp_type", "T", "s_star", "delta",
                  "method_id", "method_label", "standardizer", "hset_name", "scale_weight_name")
  pieces <- split(ok, interaction(ok[group_vars], drop = TRUE, lex.order = TRUE))
  rows <- lapply(pieces, function(z) {
    key <- z[1L, group_vars, drop = FALSE]
    key$raw_power <- mean(as.numeric(z$reject), na.rm = TRUE)
    key$sap <- mean(z$size_adjusted_reject, na.rm = TRUE)
    key$add <- mean(z$delay, na.rm = TRUE)
    key$n_valid <- sum(!is.na(z$reject))
    key$n_detected <- sum(!is.na(z$delay))
    key
  })
  do.call(rbind, rows)
}

null_par <- list(
  scenario = "null",
  dgp_type = "fMA1",
  m = 500L,
  T_value = 2,
  delta = 0,
  s_star = NA_integer_,
  basis_k = localized_change_basis_index_default(),
  alpha_levels = alpha_levels,
  q_cap = q_cap,
  fve_threshold = 0.95,
  fixed_q = NA_integer_,
  n_grid = 301L,
  ridge = 1e-10,
  range_floor = 1e-8,
  hac_bandwidth = NULL,
  page_length_grid_size = 40L,
  weighted_length_grid_size = 40L,
  finite_eval_grid_size = 250L,
  exact_page_scan = FALSE,
  exact_weighted_scan = FALSE,
  multiscale_h_sets = multiscale_h_sets
)
null_summary_file <- file.path(out_dir, "null_summary.csv")
if (file.exists(null_summary_file) && !force_null) {
  message("Reusing targeted null simulation summary: ", null_summary_file)
  null_summary <- read.csv(null_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  message("Running targeted null simulation ...")
  null_raw <- run_combo(null_par, nsim = nsim, ncores = ncores, batch_size = batch_size, cv_alt = cv_alt, alt_catalog = alt_catalog, root = ROOT)
  null_summary <- summarize_null(null_raw)
  write_csv_atomic(null_summary, null_summary_file, row.names = FALSE)
}

cv_emp <- null_summary[, c("method_id", "empirical_cv")]

combo_summaries <- list()
idx <- 1L

for (scenario_name in scenario_vec) {
  delta_vals <- delta_vec_arg %||% power_delta_grid_default(scenario_name, "fMA1", include_zero = FALSE)
  for (s_star in s_star_vec) {
    for (delta in delta_vals) {
      combo_file <- file.path(
        combo_dir,
        sprintf(
          "combo_%s_sstar%s_delta%s.csv",
          scenario_name,
          s_star,
          sanitize_tag(format(delta, trim = TRUE, scientific = FALSE))
        )
      )
      if (file.exists(combo_file)) {
        message("Reusing targeted power simulation: ", basename(combo_file))
        combo_summaries[[idx]] <- read.csv(combo_file, stringsAsFactors = FALSE, check.names = FALSE)
        idx <- idx + 1L
        next
      }
      sim_par <- list(
        scenario = scenario_name,
        dgp_type = "fMA1",
        m = 500L,
        T_value = 2,
        delta = delta,
        s_star = s_star,
        basis_k = localized_change_basis_index_default(),
        alpha_levels = alpha_levels,
        q_cap = q_cap,
        fve_threshold = 0.95,
        fixed_q = NA_integer_,
        n_grid = 301L,
        ridge = 1e-10,
        range_floor = 1e-8,
        hac_bandwidth = NULL,
        page_length_grid_size = 40L,
        weighted_length_grid_size = 40L,
        finite_eval_grid_size = 250L,
        exact_page_scan = FALSE,
        exact_weighted_scan = FALSE,
        multiscale_h_sets = multiscale_h_sets
      )
      message(sprintf("Running targeted power simulation: scenario=%s, s*= %d, delta=%s", scenario_name, s_star, format(delta, trim = TRUE)))
      raw_df <- run_combo(sim_par, nsim = nsim, ncores = ncores, batch_size = batch_size, cv_alt = cv_alt, alt_catalog = alt_catalog, root = ROOT)
      combo_summaries[[idx]] <- summarize_power_combo(raw_df, cv_emp)
      write_csv_atomic(combo_summaries[[idx]], combo_file, row.names = FALSE)
      idx <- idx + 1L
    }
  }
}

power_summary <- do.call(rbind, combo_summaries)
write_csv_atomic(power_summary, file.path(out_dir, "power_summary_by_combo.csv"), row.names = FALSE)

overall_summary <- aggregate(
  cbind(raw_power, sap, add) ~ scenario + standardizer + hset_name + scale_weight_name,
  data = power_summary,
  FUN = function(x) mean(x, na.rm = TRUE)
)
overall_summary <- overall_summary[order(overall_summary$scenario, overall_summary$standardizer, overall_summary$sap, decreasing = TRUE), ]
write_csv_atomic(overall_summary, file.path(out_dir, "overall_summary.csv"), row.names = FALSE)

best_summary <- do.call(rbind, lapply(split(overall_summary, list(overall_summary$scenario, overall_summary$standardizer), drop = TRUE), function(df) df[which.max(df$sap), , drop = FALSE]))
best_summary <- best_summary[order(best_summary$scenario, best_summary$standardizer), ]
write_csv_atomic(best_summary, file.path(out_dir, "best_by_scenario_and_standardizer.csv"), row.names = FALSE)

message("Targeted multiscale tuning diagnostic complete.")
message("Results written to: ", normalize_path2(out_dir, mustWork = FALSE))
