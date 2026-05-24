
# ==============================================================
# openend_alt_detector_critical_values.R
# ============================================================== 

g_gamma_sq <- function(s, gamma) (1 + s)^2 * (s / (1 + s))^(2 * gamma)

build_openend_alt_meta <- function(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names) {
  build_finite_alt_meta(q_grid = q_grid, gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec,
                        weighted_cusum_omega_names = weighted_cusum_omega_names,
                        multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names)
}

simulate_one_openend_alt_rep <- function(q_grid, n_train_grid = 1500L, n_open_grid = 2000L,
                                         gamma_vec = c(0, 0.15), mosum_h_vec = c(0.10, 0.20),
                                         weighted_cusum_omega_names = c("InvSqrt"),
                                         multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                         multiscale_scale_names = c("Equal"), open_pair_lag_grid_size = 40L,
                                         open_eval_grid_size = 250L, ridge = 1e-10, range_floor = 1e-8,
                                         exact_page_scan = FALSE, exact_weighted_scan = FALSE) {
  q_grid <- sort(unique(as.integer(q_grid))); q_max <- max(q_grid)
  r <- (1:n_train_grid) / n_train_grid
  dB <- matrix(rnorm(n_train_grid * q_max), nrow = n_train_grid, ncol = q_max) / sqrt(n_train_grid)
  B <- apply(dB, 2, cumsum)
  B1 <- B[n_train_grid, ]
  B0 <- B - outer(r, B1)
  V <- crossprod(B0) / n_train_grid
  ranges <- apply(B0, 2, max) - apply(B0, 2, min)
  ranges[ranges < range_floor] <- range_floor
  inv2 <- 1 / (ranges^2)

  x <- (1:n_open_grid) / (n_open_grid + 1)
  x_full <- c(0, x)
  dx <- diff(x_full)
  dG <- matrix(rnorm(n_open_grid * q_max), nrow = n_open_grid, ncol = q_max)
  dG <- dG * matrix(rep(sqrt(dx), q_max), nrow = n_open_grid, ncol = q_max)
  G <- apply(dG, 2, cumsum)
  U0_full <- rbind(rep(0, q_max), sweep(G, 1, 1 - x, FUN = "/"))
  s <- x / (1 - x)
  s_full <- c(0, s)
  eval_idx <- build_eval_index_grid(n_open_grid, open_eval_grid_size)

  page_lag_grid <- full_or_geometric_lag_grid(n_open_grid, open_pair_lag_grid_size, exact = exact_page_scan)
  weighted_lag_grid <- full_or_geometric_lag_grid(n_open_grid, open_pair_lag_grid_size, exact = exact_weighted_scan)
  store <- empty_value_store(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names)

  norms_full <- compute_metric_norms(U0_full[eval_idx + 1L, , drop = FALSE], V = V, inv2 = inv2, ridge = ridge)
  s_eval <- s[eval_idx]
  for (std in names(norms_full)) {
    M <- norms_full[[std]][, q_grid, drop = FALSE]
    for (g in gamma_vec) store$full_vals[[std]][[as.character(g)]] <- apply_scale_and_max(M, 1 / g_gamma_sq(s_eval, g))
  }

  for (L in page_lag_grid) {
    hi <- eval_idx[eval_idx >= L]
    if (length(hi) == 0L) next
    lo <- hi - L
    D <- U0_full[hi + 1L, , drop = FALSE] - U0_full[lo + 1L, , drop = FALSE]
    lens <- s_full[hi + 1L] - s_full[lo + 1L]
    norms <- compute_metric_norms(D, V = V, inv2 = inv2, ridge = ridge)
    for (std in names(norms)) {
      M <- norms[[std]][, q_grid, drop = FALSE]
      for (g in gamma_vec) {
        glab <- as.character(g)
        store$page_vals[[std]][[glab]] <- pmax(store$page_vals[[std]][[glab]], apply_scale_and_max(M, 1 / g_gamma_sq(lens, g)))
      }
    }
  }

  for (L in weighted_lag_grid) {
    hi <- eval_idx[eval_idx >= L]
    if (length(hi) == 0L) next
    lo <- hi - L
    D <- U0_full[hi + 1L, , drop = FALSE] - U0_full[lo + 1L, , drop = FALSE]
    lens <- s_full[hi + 1L] - s_full[lo + 1L]
    norms <- compute_metric_norms(D, V = V, inv2 = inv2, ridge = ridge)
    for (std in names(norms)) {
      M <- norms[[std]][, q_grid, drop = FALSE]
      for (g in gamma_vec) {
        glab <- as.character(g)
        base_sc <- 1 / g_gamma_sq(lens, g)
        for (om in weighted_cusum_omega_names) {
          sc <- base_sc * make_weighted_cusum_omega(lens, om)
          store$weighted_vals[[std]][[glab]][[om]] <- pmax(store$weighted_vals[[std]][[glab]][[om]], apply_scale_and_max(M, sc))
        }
      }
    }
  }

  for (h in mosum_h_vec) {
    valid_hi <- eval_idx[s_full[eval_idx + 1L] >= 2 * h]
    if (length(valid_hi) == 0L) next
    current <- U0_full[valid_hi + 1L, , drop = FALSE]
    targ1 <- s_full[valid_hi + 1L] - h
    targ2 <- s_full[valid_hi + 1L] - 2 * h
    interp1 <- interp_matrix(s_full, U0_full, targ1)
    interp2 <- interp_matrix(s_full, U0_full, targ2)
    D <- current - 2 * interp1 + interp2
    norms <- compute_metric_norms(D, V = V, inv2 = inv2, ridge = ridge)
    for (std in names(norms)) {
      M <- norms[[std]][, q_grid, drop = FALSE] / h
      hlab <- format(h, trim = TRUE)
      cur <- col_maxs_base(M)
      store$mosum_vals[[std]][[hlab]] <- pmax(store$mosum_vals[[std]][[hlab]], cur)
      for (hs in names(multiscale_h_sets)) {
        if (any(abs(multiscale_h_sets[[hs]] - h) < 1e-10)) {
          for (an in multiscale_scale_names) {
            store$multi_vals[[std]][[hs]][[an]] <- pmax(store$multi_vals[[std]][[hs]][[an]], cur * make_multiscale_weight(h, an))
          }
        }
      }
    }
  }

  values <- numeric(0L)
  for (q in q_grid) {
    qpos <- match(q, q_grid)
    for (std in c("HAC", "SSMS", "RSMS")) {
      for (g in gamma_vec) {
        glab <- as.character(g)
        values <- c(values, store$full_vals[[std]][[glab]][qpos])
        values <- c(values, store$page_vals[[std]][[glab]][qpos])
      }
      for (h in mosum_h_vec) values <- c(values, store$mosum_vals[[std]][[format(h, trim = TRUE)]][qpos])
      for (g in gamma_vec) for (om in weighted_cusum_omega_names) values <- c(values, store$weighted_vals[[std]][[as.character(g)]][[om]][qpos])
      for (hs in names(multiscale_h_sets)) for (an in multiscale_scale_names) values <- c(values, store$multi_vals[[std]][[hs]][[an]][qpos])
    }
  }
  values
}

simulate_openend_alt_critical_values <- function(q_grid, gamma_vec = c(0, 0.15), mosum_h_vec = c(0.10, 0.20),
                                                 weighted_cusum_omega_names = c("InvSqrt"),
                                                 multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                                 multiscale_scale_names = c("Equal"), alpha_levels = c(0.10, 0.05, 0.01),
                                                 nrep = 5000L, n_train_grid = 1500L, n_open_grid = 2000L,
                                                 open_pair_lag_grid_size = 40L, open_eval_grid_size = 250L,
                                                 ridge = 1e-10, range_floor = 1e-8, ncores = 1L, seed = 13579L,
                                                 exact_page_scan = FALSE, exact_weighted_scan = FALSE, batch_size = 100L,
                                                 verbose = TRUE) {
  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  mosum_h_vec <- sort(unique(as.numeric(mosum_h_vec)))
  weighted_cusum_omega_names <- vapply(weighted_cusum_omega_names, normalize_omega_name, character(1L))
  multiscale_scale_names <- vapply(multiscale_scale_names, normalize_scale_weight_name, character(1L))
  alpha_levels <- as.numeric(alpha_levels)
  meta <- build_openend_alt_meta(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names)
  seeds <- make_seed_stream(nrep = nrep, seed = seed)
  core_plan <- resolve_safe_ncores(requested = ncores)
  ncores_use <- min(core_plan$used, batch_size)
  batches <- split(seeds, ceiling(seq_along(seeds) / batch_size))
  res_parts <- vector("list", length(batches))

  if (verbose) message(sprintf("Starting open-end alternative-detector simulation: T=Inf | q-grid=%s | nrep=%d | train-grid=%d | open-grid=%d | workers(requested=%d, used=%d) | batches=%d", format_q_grid(q_grid), nrep, n_train_grid, n_open_grid, core_plan$requested, ncores_use, length(batches)))

  for (b in seq_along(batches)) {
    seed_batch <- batches[[b]]
    if (ncores_use <= 1L || length(seed_batch) == 1L) {
      res_batch <- lapply(seed_batch, function(sd) {
        set.seed(sd)
        simulate_one_openend_alt_rep(q_grid = q_grid, n_train_grid = n_train_grid, n_open_grid = n_open_grid,
          gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec, weighted_cusum_omega_names = weighted_cusum_omega_names,
          multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names,
          open_pair_lag_grid_size = open_pair_lag_grid_size, open_eval_grid_size = open_eval_grid_size,
          ridge = ridge, range_floor = range_floor, exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan)
      })
    } else {
      cl <- make_psock_cluster(min(ncores_use, length(seed_batch)))
      on.exit(stop_psock_cluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c("simulate_one_openend_alt_rep", "compute_metric_norms", "empty_value_store",
        "row_cumsum_weighted_squares", "leading_mahalanobis_sequence", "apply_scale_and_max", "col_maxs_base",
        "build_eval_index_grid", "build_linear_index_grid", "build_geometric_index_grid", "full_or_geometric_lag_grid", "interp_matrix",
        "g_gamma_sq", "make_weighted_cusum_omega", "make_multiscale_weight", "normalize_omega_name",
        "normalize_scale_weight_name", "q_grid", "n_train_grid", "n_open_grid", "gamma_vec", "mosum_h_vec",
        "weighted_cusum_omega_names", "multiscale_h_sets", "multiscale_scale_names", "open_pair_lag_grid_size",
        "open_eval_grid_size", "ridge", "range_floor", "exact_page_scan", "exact_weighted_scan"), envir = environment())
      res_batch <- parallel::parLapplyLB(cl, seed_batch, function(sd) {
        set.seed(sd)
        simulate_one_openend_alt_rep(q_grid = q_grid, n_train_grid = n_train_grid, n_open_grid = n_open_grid,
          gamma_vec = gamma_vec, mosum_h_vec = mosum_h_vec, weighted_cusum_omega_names = weighted_cusum_omega_names,
          multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names,
          open_pair_lag_grid_size = open_pair_lag_grid_size, open_eval_grid_size = open_eval_grid_size,
          ridge = ridge, range_floor = range_floor, exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan)
      })
      stop_psock_cluster(cl)
      on.exit(NULL, add = FALSE)
    }
    res_parts[[b]] <- do.call(rbind, res_batch)
    if (verbose) message(sprintf("  T=Inf batch %d/%d finished (%d/%d reps)", b, length(batches), min(b * batch_size, nrep), nrep))
  }

  res_mat <- do.call(rbind, res_parts)
  source_tag <- sprintf("simulated_alt_detectors_openend_q_%d_%d_nrep_%d.csv", min(q_grid), max(q_grid), nrep)
  rows <- list(); idx <- 1L
  for (a in alpha_levels) {
    for (j in seq_len(nrow(meta))) {
      rowj <- meta[j, , drop = FALSE]
      rows[[idx]] <- data.frame(standardizer = rowj$standardizer, detector = rowj$detector, T = "Inf", gamma = rowj$gamma,
        q = rowj$q, alpha = a, critical_value = quantile_upper(res_mat[, j], a), source_file = source_tag,
        bandwidth_h = rowj$bandwidth_h, omega_name = rowj$omega_name, hset_name = rowj$hset_name,
        scale_weight_name = rowj$scale_weight_name, exact_scan = FALSE, stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  if (verbose) message("Finished open-end alternative-detector simulation for T=Inf.")
  list(meta = meta, rows = if (length(rows)) do.call(rbind, rows) else empty_alt_table(), info = list(T = Inf, q_grid = q_grid, nrep = nrep, n_train_grid = n_train_grid, n_open_grid = n_open_grid, workers_used = ncores_use))
}
