# ==============================================================
# utils.R -- helpers for simulation, reporting, and tables
# ============================================================== 

is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

safe_solve <- function(A, ridge = 1e-10) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  q <- nrow(A)
  A2 <- A + ridge * diag(q)

  tryCatch(
    chol2inv(chol(A2)),
    error = function(e) solve(A2)
  )
}

get_connection_limit <- function() {
  opt <- getOption("connections")
  if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt >= 16) {
    return(as.integer(opt))
  }

  envv <- suppressWarnings(as.integer(Sys.getenv("R_CONNECTIONS_LIMIT", unset = "")))
  if (!is.na(envv) && envv >= 16L) {
    return(envv)
  }

  128L
}

current_open_connections <- function() {
  out <- tryCatch(showConnections(all = TRUE), error = function(e) NULL)
  if (is.null(out)) return(0L)
  as.integer(nrow(out))
}

resolve_safe_ncores <- function(requested = NULL,
                                reserve_connections = 8L,
                                hard_cap_windows = 32L,
                                hard_cap_other = 64L) {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1L) detected <- 1L
  detected <- as.integer(detected)

  requested_default <- max(1L, detected - 1L)
  requested_value <- if (is.null(requested) || length(requested) == 0L || is.na(requested)) {
    requested_default
  } else {
    max(1L, as.integer(requested[1L]))
  }

  connection_limit <- get_connection_limit()
  open_connections <- current_open_connections()
  connection_budget <- max(1L, connection_limit - open_connections - as.integer(reserve_connections))
  os_cap <- if (is_windows()) as.integer(hard_cap_windows) else as.integer(hard_cap_other)

  used <- min(requested_value, detected, connection_budget, os_cap)
  used <- max(1L, as.integer(used))

  reasons <- character(0L)
  if (requested_value > detected) reasons <- c(reasons, "requested exceeds detected cores")
  if (requested_value > connection_budget) reasons <- c(reasons, "requested exceeds available connection budget")
  if (requested_value > os_cap) reasons <- c(reasons, if (is_windows()) "requested exceeds Windows safety cap" else "requested exceeds safety cap")
  if (length(reasons) == 0L) reasons <- "no cap applied"

  list(
    requested = as.integer(requested_value),
    detected = as.integer(detected),
    used = as.integer(used),
    connection_limit = as.integer(connection_limit),
    open_connections = as.integer(open_connections),
    reserve_connections = as.integer(reserve_connections),
    connection_budget = as.integer(connection_budget),
    os_cap = as.integer(os_cap),
    reasons = reasons
  )
}

make_psock_cluster <- function(ncores) {
  ncores <- as.integer(max(1L, ncores))
  if (ncores <= 1L) return(NULL)
  parallel::makeCluster(ncores, type = "PSOCK", outfile = "")
}

stop_psock_cluster <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  invisible(NULL)
}

split_into_chunks <- function(x, nchunks) {
  n <- length(x)
  if (n == 0L) return(list())
  nchunks <- max(1L, min(as.integer(nchunks), n))
  split(x, cut(seq_len(n), breaks = nchunks, labels = FALSE))
}

canonical_T_scalar <- function(T) {
  if (length(T) != 1L) stop("T must be scalar.")

  if (is.character(T)) {
    x <- trimws(T)
    if (tolower(x) %in% c("inf", "infinity")) return("Inf")
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)) T <- num else return(x)
  }

  if (is.infinite(T)) return("Inf")
  if (!is.numeric(T) || is.na(T)) stop("Invalid T value.")

  if (abs(T - round(T)) < 1e-10) {
    return(as.character(as.integer(round(T))))
  }
  as.character(T)
}

canonical_T_vec <- function(T) {
  vapply(T, canonical_T_scalar, character(1L))
}

row_match_numeric <- function(x, target, tol = 1e-10) {
  if (is.na(target)) return(rep(TRUE, length(x)))
  if (all(is.na(x))) return(rep(TRUE, length(x)))
  (!is.na(x)) & (abs(x - target) < tol)
}

trapz_equal_grid <- function(y) {
  n <- length(y)
  if (n < 2L) return(0)
  dx <- 1 / (n - 1)
  if (n == 2L) return(dx * (0.5 * y[1L] + 0.5 * y[2L]))
  dx * (0.5 * y[1L] + sum(y[2L:(n - 1L)]) + 0.5 * y[n])
}

quantile_upper <- function(x, alpha) {
  as.numeric(stats::quantile(x, probs = 1 - alpha, names = FALSE, type = 7))
}

make_seed_stream <- function(nrep, seed = 123456789L) {
  max_int <- .Machine$integer.max
  seeds <- (as.integer(seed) + seq_len(nrep)) %% max_int
  seeds[seeds == 0L] <- 1L
  seeds
}

empty_base_table <- function() {
  data.frame(
    stat = character(0L),
    type = character(0L),
    T = character(0L),
    gamma = numeric(0L),
    q = integer(0L),
    alpha = numeric(0L),
    critical_value = numeric(0L),
    source_file = character(0L),
    stringsAsFactors = FALSE
  )
}

empty_weight_table <- function() {
  data.frame(
    stat = character(0L),
    type = character(0L),
    T = character(0L),
    gamma = numeric(0L),
    q = integer(0L),
    alpha = numeric(0L),
    critical_value = numeric(0L),
    source_file = character(0L),
    weight_name = character(0L),
    stringsAsFactors = FALSE
  )
}

format_q_grid <- function(q_grid) {
  q_grid <- sort(unique(as.integer(q_grid)))
  if (length(q_grid) == 0L) return("<empty>")
  if (length(q_grid) <= 10L) return(paste(q_grid, collapse = ", "))
  sprintf("%d:%d (%d values)", min(q_grid), max(q_grid), length(q_grid))
}

format_num_vec <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

report_simulation_design <- function(root,
                                     output_dir,
                                     q_grid,
                                     T_grid,
                                     gamma_vec,
                                     weight_names,
                                     alpha_levels,
                                     nrep_finite,
                                     nrep_openend,
                                     grid_monitor,
                                     n_train_grid_open,
                                     n_open_grid_open,
                                     core_plan,
                                     overwrite_existing = TRUE) {
  lines <- c(
    strrep("=", 78L),
    "Streaming-curve critical-value simulation",
    "Fresh run: all requested critical values will be simulated from scratch.",
    sprintf("Project root: %s", normalize_path2(root, mustWork = FALSE)),
    sprintf("Output directory: %s", normalize_path2(output_dir, mustWork = FALSE)),
    "Statistics:",
    "  - Shao's KS (SSMS, KS)",
    "  - Adjusted-range based KS (RSMS, KS)",
    "  - HAC-based KS (HAC, KS)",
    "  - Shao's weighted CvM (SSMS, CvM)",
    "  - Adjusted-range based weighted CvM (RSMS, CvM)",
    "  - HAC-based weighted CvM (HAC, CvM)",
    sprintf("q-grid: %s", format_q_grid(q_grid)),
    sprintf("T-grid: %s", paste(canonical_T_vec(T_grid), collapse = ", ")),
    sprintf("gamma values for KS: %s", format_num_vec(gamma_vec)),
    sprintf("CvM weights: %s", paste(vapply(weight_names, normalize_weight_name, character(1L)), collapse = ", ")),
    sprintf("alpha levels: %s", format_num_vec(alpha_levels)),
    sprintf("Finite-horizon replications: %d", as.integer(nrep_finite)),
    sprintf("Open-end replications: %d", as.integer(nrep_openend)),
    sprintf("Finite-horizon grid_monitor: %d", as.integer(grid_monitor)),
    sprintf("Open-end training grid size: %d", as.integer(n_train_grid_open)),
    sprintf("Open-end monitoring grid size: %d", as.integer(n_open_grid_open)),
    sprintf("Parallel workers requested: %d", as.integer(core_plan$requested)),
    sprintf("Detected logical cores: %d", as.integer(core_plan$detected)),
    sprintf("Open connections at start: %d", as.integer(core_plan$open_connections)),
    sprintf("Connection limit used: %d", as.integer(core_plan$connection_limit)),
    sprintf("Connection reserve: %d", as.integer(core_plan$reserve_connections)),
    sprintf("Connection-budget cap: %d", as.integer(core_plan$connection_budget)),
    sprintf("OS safety cap: %d", as.integer(core_plan$os_cap)),
    sprintf("Parallel workers to be used: %d", as.integer(core_plan$used)),
    sprintf("Parallel core rule: %s", paste(core_plan$reasons, collapse = "; ")),
    sprintf("Existing output files will %s.", if (overwrite_existing) "be backed up and overwritten" else "be preserved"),
    "Implied finite-horizon (m, k_max) pairs:")

  finite_T <- T_grid[is.finite(T_grid)]
  if (length(finite_T) == 0L) {
    lines <- c(lines, "  - none")
  } else {
    for (T in finite_T) {
      m <- as.integer(round(grid_monitor / T))
      k_max <- as.integer(round(m * T))
      lines <- c(lines, sprintf("  - T=%s: m=%d, k_max=%d", canonical_T_scalar(T), m, k_max))
    }
  }

  lines <- c(lines, strrep("=", 78L))
  message(paste(lines, collapse = "\n"))
  invisible(lines)
}

write_run_summary <- function(lines, output_dir, filename = "run_summary.txt") {
  ensure_dir(output_dir)
  path <- file.path(output_dir, filename)
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}
