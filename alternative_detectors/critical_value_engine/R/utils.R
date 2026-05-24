
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
  tryCatch(chol2inv(chol(A2)), error = function(e) solve(A2))
}

safe_chol <- function(A, ridge = 1e-10) {
  q <- nrow(A)
  A2 <- A + ridge * diag(q)
  tryCatch(chol(A2), error = function(e) chol(A2 + 10 * ridge * diag(q)))
}

get_connection_limit <- function() {
  opt <- getOption("connections")
  if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt >= 16) return(as.integer(opt))
  envv <- suppressWarnings(as.integer(Sys.getenv("R_CONNECTIONS_LIMIT", unset = "")))
  if (!is.na(envv) && envv >= 16L) return(envv)
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
  requested_value <- if (is.null(requested) || length(requested) == 0L || is.na(requested)) requested_default else max(1L, as.integer(requested[1L]))
  connection_limit <- get_connection_limit()
  open_connections <- current_open_connections()
  connection_budget <- max(1L, connection_limit - open_connections - as.integer(reserve_connections))
  os_cap <- if (is_windows()) as.integer(hard_cap_windows) else as.integer(hard_cap_other)
  used <- max(1L, min(requested_value, detected, connection_budget, os_cap))
  reasons <- character(0L)
  if (requested_value > detected) reasons <- c(reasons, "requested exceeds detected cores")
  if (requested_value > connection_budget) reasons <- c(reasons, "requested exceeds available connection budget")
  if (requested_value > os_cap) reasons <- c(reasons, if (is_windows()) "requested exceeds Windows safety cap" else "requested exceeds safety cap")
  if (length(reasons) == 0L) reasons <- "no cap applied"
  list(requested = as.integer(requested_value), detected = as.integer(detected), used = as.integer(used),
       connection_limit = as.integer(connection_limit), open_connections = as.integer(open_connections),
       reserve_connections = as.integer(reserve_connections), connection_budget = as.integer(connection_budget),
       os_cap = as.integer(os_cap), reasons = reasons)
}

make_psock_cluster <- function(ncores) {
  ncores <- as.integer(max(1L, ncores))
  if (ncores <= 1L) return(NULL)
  parallel::makeCluster(ncores, type = "PSOCK", outfile = "")
}

stop_psock_cluster <- function(cl) {
  if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  invisible(NULL)
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
  if (abs(T - round(T)) < 1e-10) return(as.character(as.integer(round(T))))
  as.character(T)
}

canonical_T_vec <- function(T) vapply(T, canonical_T_scalar, character(1L))

row_match_numeric <- function(x, target, tol = 1e-10) {
  if (is.na(target)) return(rep(TRUE, length(x)))
  if (all(is.na(x))) return(rep(TRUE, length(x)))
  (!is.na(x)) & (abs(x - target) < tol)
}

quantile_upper <- function(x, alpha) as.numeric(stats::quantile(x, probs = 1 - alpha, names = FALSE, type = 7))

make_seed_stream <- function(nrep, seed = 123456789L) {
  max_int <- .Machine$integer.max
  seeds <- (as.integer(seed) + seq_len(nrep)) %% max_int
  seeds[seeds == 0L] <- 1L
  seeds
}

format_q_grid <- function(q_grid) {
  q_grid <- sort(unique(as.integer(q_grid)))
  if (length(q_grid) == 0L) return("<empty>")
  if (length(q_grid) <= 10L) return(paste(q_grid, collapse = ", "))
  sprintf("%d:%d (%d values)", min(q_grid), max(q_grid), length(q_grid))
}

format_num_vec <- function(x) paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
format_h_vec <- function(x) paste(format(signif(x, 6), trim = TRUE, scientific = FALSE), collapse = ", ")

build_geometric_index_grid <- function(n_max, n_points = 250L) {
  n_max <- as.integer(max(1L, n_max))
  n_points <- as.integer(max(1L, min(n_points, n_max)))
  if (n_points >= n_max) return(seq_len(n_max))
  out <- unique(c(1L, as.integer(round(exp(seq(log(1), log(n_max), length.out = n_points)))), n_max))
  out[out >= 1L & out <= n_max]
}

build_linear_index_grid <- function(n_max, n_points = 250L) {
  n_max <- as.integer(max(1L, n_max))
  n_points <- as.integer(max(1L, min(n_points, n_max)))
  if (n_points >= n_max) return(seq_len(n_max))
  out <- unique(c(1L, as.integer(round(seq(1, n_max, length.out = n_points))), n_max))
  out[out >= 1L & out <= n_max]
}

build_eval_index_grid <- function(n_max, n_points = 250L) {
  n_max <- as.integer(max(1L, n_max))
  n_points <- as.integer(max(1L, min(n_points, n_max)))
  if (n_points >= n_max) return(seq_len(n_max))
  out <- unique(sort(c(
    1L,
    seq_len(min(25L, n_max)),
    build_linear_index_grid(n_max, n_points),
    build_geometric_index_grid(n_max, n_points),
    n_max
  )))
  out[out >= 1L & out <= n_max]
}

full_or_geometric_lag_grid <- function(n_max, grid_size = 250L, exact = FALSE) {
  n_max <- as.integer(max(1L, n_max))
  if (isTRUE(exact) || isTRUE(grid_size >= n_max)) return(seq_len(n_max))
  build_geometric_index_grid(n_max = n_max, n_points = grid_size)
}

interp_matrix <- function(x, Y, xout) {
  stopifnot(length(x) == nrow(Y))
  x <- as.numeric(x); xout <- as.numeric(xout)
  out <- matrix(NA_real_, nrow = length(xout), ncol = ncol(Y))
  for (j in seq_len(ncol(Y))) {
    out[, j] <- approx(x = x, y = Y[, j], xout = xout, method = "linear", ties = mean, rule = 2)$y
  }
  out
}

row_cumsum_weighted_squares <- function(X, w = NULL) {
  X2 <- X^2
  if (!is.null(w)) X2 <- sweep(X2, 2L, w, `*`)
  out <- X2
  q <- ncol(X2)
  if (q >= 2L) {
    for (j in 2:q) out[, j] <- out[, j - 1L] + X2[, j]
  }
  out
}

leading_mahalanobis_sequence <- function(X, V, ridge = 1e-10) {
  n <- nrow(X); qmax <- ncol(X)
  out <- matrix(0, n, qmax)
  A11 <- as.numeric(V[1L, 1L])
  if (!is.finite(A11) || A11 <= ridge) A11 <- ridge
  Ainv <- matrix(1 / A11, nrow = 1L, ncol = 1L)
  out[, 1L] <- X[, 1L]^2 * Ainv[1L, 1L]
  if (qmax == 1L) return(out)
  for (q in 2:qmax) {
    b <- V[1:(q - 1L), q, drop = FALSE]
    c0 <- as.numeric(V[q, q])
    Ainv_b <- Ainv %*% b
    schur <- c0 - as.numeric(crossprod(b, Ainv_b))
    if (!is.finite(schur) || schur <= ridge) schur <- ridge
    resid <- X[, q] - X[, 1:(q - 1L), drop = FALSE] %*% Ainv_b
    out[, q] <- out[, q - 1L] + as.numeric(resid)^2 / schur
    top_left <- Ainv + (Ainv_b %*% t(Ainv_b)) / schur
    col_new <- -Ainv_b / schur
    Ainv <- rbind(cbind(top_left, col_new), cbind(t(col_new), 1 / schur))
  }
  out
}

col_maxs_base <- function(M) {
  if (is.null(dim(M))) return(as.numeric(M))
  if (nrow(M) == 1L) return(as.numeric(M[1L, ]))
  apply(M, 2L, max)
}

apply_scale_and_max <- function(M, scale_vec) {
  col_maxs_base(sweep(M, 1L, scale_vec, `*`))
}

empty_alt_table <- function() {
  data.frame(
    standardizer = character(0L), detector = character(0L), T = character(0L), gamma = numeric(0L),
    q = integer(0L), alpha = numeric(0L), critical_value = numeric(0L), source_file = character(0L),
    bandwidth_h = numeric(0L), omega_name = character(0L), hset_name = character(0L),
    scale_weight_name = character(0L), exact_scan = logical(0L), stringsAsFactors = FALSE
  )
}

report_simulation_design_alt <- function(root, output_dir, q_grid, T_grid, gamma_vec, alpha_levels,
                                         nrep_finite, nrep_openend, grid_monitor, n_train_grid_open,
                                         n_open_grid_open, core_plan, mosum_h_vec, weighted_cusum_omega_names,
                                         multiscale_h_sets, multiscale_scale_names, page_length_grid_size,
                                         weighted_length_grid_size, finite_eval_grid_size, open_eval_grid_size,
                                         batch_size, overwrite_existing = TRUE) {
  hset_lines <- if (length(multiscale_h_sets) == 0L) "  - none" else vapply(names(multiscale_h_sets), function(nm) sprintf("  - %s: %s", nm, format_h_vec(multiscale_h_sets[[nm]])), character(1L))
  c(
    strrep("=", 78L),
    "Alternative-detector critical-value simulation",
    "Fresh run: all requested critical values will be simulated from scratch.",
    "This fast version uses thinned scan grids for benchmark detectors so the job completes in plain R.",
    sprintf("Project root: %s", normalize_path2(root, mustWork = FALSE)),
    sprintf("Output directory: %s", normalize_path2(output_dir, mustWork = FALSE)),
    "Detector families:", "  - Full-CUSUM", "  - Page-CUSUM", "  - MOSUM", "  - weighted-CUSUM", "  - multiscale MOSUM",
    "Standardizers:", "  - HAC", "  - SSMS (Shao self-normalization)", "  - RSMS (adjusted-range self-normalization)",
    sprintf("q-grid: %s", format_q_grid(q_grid)),
    sprintf("T-grid: %s", paste(canonical_T_vec(T_grid), collapse = ", ")),
    sprintf("gamma values: %s", format_num_vec(gamma_vec)),
    sprintf("alpha levels: %s", format_num_vec(alpha_levels)),
    sprintf("finite-horizon replications: %d", nrep_finite),
    sprintf("open-end replications: %d", nrep_openend),
    sprintf("finite grid_monitor: %d", grid_monitor),
    sprintf("open-end training grid: %d", n_train_grid_open),
    sprintf("open-end monitoring grid: %d", n_open_grid_open),
    sprintf("finite evaluation grid size: %d", finite_eval_grid_size),
    sprintf("open-end evaluation grid size: %d", open_eval_grid_size),
    sprintf("Page lag-grid size: %d", page_length_grid_size),
    sprintf("weighted-CUSUM lag-grid size: %d", weighted_length_grid_size),
    sprintf("MOSUM bandwidth h values: %s", format_h_vec(mosum_h_vec)),
    sprintf("weighted-CUSUM omega choices: %s", paste(weighted_cusum_omega_names, collapse = ", ")),
    "multiscale h-sets:", hset_lines,
    sprintf("multiscale scale weights: %s", paste(multiscale_scale_names, collapse = ", ")),
    sprintf("batch size: %d", batch_size),
    sprintf("workers requested: %d", core_plan$requested),
    sprintf("workers used: %d", core_plan$used),
    sprintf("detected cores: %d", core_plan$detected),
    sprintf("connection budget: %d", core_plan$connection_budget),
    sprintf("OS safety cap: %d", core_plan$os_cap),
    sprintf("fresh overwrite existing outputs: %s", if (overwrite_existing) "YES" else "NO"),
    strrep("=", 78L)
  )
}

write_run_summary <- function(lines, out_dir, filename = "run_summary.txt") {
  ensure_dir(out_dir)
  path <- file.path(out_dir, filename)
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}
