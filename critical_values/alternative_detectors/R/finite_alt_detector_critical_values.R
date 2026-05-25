
# ==============================================================
# finite_alt_detector_critical_values.R
# ============================================================== 

compute_m_from_T <- function(grid_monitor = 10000L, T) {
  if (!is.finite(T)) stop("compute_m_from_T() is for finite T only.")
  m <- grid_monitor / T
  if (abs(m - round(m)) > 1e-8) stop(sprintf("grid_monitor=%s is not divisible by T=%s.", grid_monitor, T))
  as.integer(round(m))
}

g_gamma_sq <- function(s, gamma) (1 + s)^2 * (s / (1 + s))^(2 * gamma)

build_finite_alt_meta <- function(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names) {
  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  mosum_h_vec <- sort(unique(as.numeric(mosum_h_vec)))
  weighted_cusum_omega_names <- vapply(weighted_cusum_omega_names, normalize_omega_name, character(1L))
  multiscale_scale_names <- vapply(multiscale_scale_names, normalize_scale_weight_name, character(1L))
  if (is.null(names(multiscale_h_sets)) || any(names(multiscale_h_sets) == "")) stop("multiscale_h_sets must be a named list.")
  rows <- list(); idx <- 1L; stdizers <- c("HAC", "SSMS", "RSMS")
  for (q in q_grid) {
    for (std in stdizers) {
      for (g in gamma_vec) {
        rows[[idx]] <- data.frame(standardizer = std, detector = "FullCUSUM", gamma = g, q = q, bandwidth_h = NA_real_, omega_name = "", hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
        rows[[idx]] <- data.frame(standardizer = std, detector = "PageCUSUM", gamma = g, q = q, bandwidth_h = NA_real_, omega_name = "", hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      }
      for (h in mosum_h_vec) {
        rows[[idx]] <- data.frame(standardizer = std, detector = "MOSUM", gamma = NA_real_, q = q, bandwidth_h = h, omega_name = "", hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      }
      for (g in gamma_vec) {
        for (om in weighted_cusum_omega_names) {
          rows[[idx]] <- data.frame(standardizer = std, detector = "WeightedCUSUM", gamma = g, q = q, bandwidth_h = NA_real_, omega_name = om, hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
        }
      }
      for (hs in names(multiscale_h_sets)) {
        for (an in multiscale_scale_names) {
          rows[[idx]] <- data.frame(standardizer = std, detector = "MultiscaleMOSUM", gamma = NA_real_, q = q, bandwidth_h = NA_real_, omega_name = "", hset_name = hs, scale_weight_name = an, stringsAsFactors = FALSE); idx <- idx + 1L
        }
      }
    }
  }
  do.call(rbind, rows)
}

compute_metric_norms <- function(D, V, inv2, ridge = 1e-10) {
  list(
    HAC = row_cumsum_weighted_squares(D),
    SSMS = leading_mahalanobis_sequence(D, V = V, ridge = ridge),
    RSMS = row_cumsum_weighted_squares(D, w = inv2)
  )
}

empty_value_store <- function(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names) {
  full_vals <- lapply(c("HAC", "SSMS", "RSMS"), function(z) setNames(vector("list", length(gamma_vec)), as.character(gamma_vec))); names(full_vals) <- c("HAC", "SSMS", "RSMS")
  page_vals <- lapply(c("HAC", "SSMS", "RSMS"), function(z) setNames(vector("list", length(gamma_vec)), as.character(gamma_vec))); names(page_vals) <- c("HAC", "SSMS", "RSMS")
  weighted_vals <- lapply(c("HAC", "SSMS", "RSMS"), function(z) setNames(lapply(as.character(gamma_vec), function(.) setNames(vector("list", length(weighted_cusum_omega_names)), weighted_cusum_omega_names)), as.character(gamma_vec))); names(weighted_vals) <- c("HAC", "SSMS", "RSMS")
  mosum_vals <- lapply(c("HAC", "SSMS", "RSMS"), function(z) setNames(vector("list", length(mosum_h_vec)), format(mosum_h_vec, trim = TRUE))); names(mosum_vals) <- c("HAC", "SSMS", "RSMS")
  multi_vals <- lapply(c("HAC", "SSMS", "RSMS"), function(z) { tmp <- list(); for (hs in names(multiscale_h_sets)) tmp[[hs]] <- setNames(vector("list", length(multiscale_scale_names)), multiscale_scale_names); tmp }); names(multi_vals) <- c("HAC", "SSMS", "RSMS")
  for (std in c("HAC", "SSMS", "RSMS")) {
    for (g in as.character(gamma_vec)) {
      full_vals[[std]][[g]] <- rep(-Inf, length(q_grid))
      page_vals[[std]][[g]] <- rep(-Inf, length(q_grid))
      for (om in weighted_cusum_omega_names) weighted_vals[[std]][[g]][[om]] <- rep(-Inf, length(q_grid))
    }
    for (h in format(mosum_h_vec, trim = TRUE)) mosum_vals[[std]][[h]] <- rep(-Inf, length(q_grid))
    for (hs in names(multiscale_h_sets)) for (an in multiscale_scale_names) multi_vals[[std]][[hs]][[an]] <- rep(-Inf, length(q_grid))
  }
  list(full_vals = full_vals, page_vals = page_vals, weighted_vals = weighted_vals, mosum_vals = mosum_vals, multi_vals = multi_vals)
}

simulate_one_finite_alt_rep <- function(q_grid, m, k_max, gamma_vec, mosum_h_vec, weighted_cusum_omega_names,
                                        multiscale_h_sets, multiscale_scale_names, page_length_grid_size,
                                        weighted_length_grid_size, finite_eval_grid_size = 250L,
                                        ridge = 1e-10, range_floor = 1e-8,
                                        exact_page_scan = FALSE, exact_weighted_scan = FALSE) {
  q_grid <- sort(unique(as.integer(q_grid)))
  q_max <- max(q_grid)
  total_steps <- m + k_max
  dW <- matrix(rnorm(total_steps * q_max), nrow = total_steps, ncol = q_max) / sqrt(m)
  W <- apply(dW, 2, cumsum)

  r <- (1:m) / m
  s <- (1:k_max) / m
  s_full <- c(0, s)
  B1 <- W[m, ]
  B0 <- W[1:m, , drop = FALSE] - outer(r, B1)
  V <- crossprod(B0) / m
  ranges <- apply(B0, 2, max) - apply(B0, 2, min)
  ranges[ranges < range_floor] <- range_floor
  inv2 <- 1 / (ranges^2)

  U_full <- W[(m + 1):(m + k_max), , drop = FALSE] - outer(1 + s, B1)
  U0_full <- rbind(rep(0, q_max), U_full)
  eval_idx <- build_eval_index_grid(k_max, finite_eval_grid_size)

  page_lag_grid <- full_or_geometric_lag_grid(k_max, page_length_grid_size, exact = exact_page_scan)
  weighted_lag_grid <- full_or_geometric_lag_grid(k_max, weighted_length_grid_size, exact = exact_weighted_scan)
  store <- empty_value_store(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names)

  # Full-CUSUM on thinned endpoint grid
  norms_full <- compute_metric_norms(U0_full[eval_idx + 1L, , drop = FALSE], V = V, inv2 = inv2, ridge = ridge)
  s_eval <- s[eval_idx]
  for (std in names(norms_full)) {
    M <- norms_full[[std]][, q_grid, drop = FALSE]
    for (g in gamma_vec) {
      store$full_vals[[std]][[as.character(g)]] <- apply_scale_and_max(M, 1 / g_gamma_sq(s_eval, g))
    }
  }

  # Page-CUSUM
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
        sc <- 1 / g_gamma_sq(lens, g)
        store$page_vals[[std]][[glab]] <- pmax(store$page_vals[[std]][[glab]], apply_scale_and_max(M, sc))
      }
    }
  }

  # Weighted-CUSUM
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

  # MOSUM / multiscale MOSUM
  for (h in mosum_h_vec) {
    G <- as.integer(floor(m * h))
    if (G < 1L || 2L * G > k_max) next
    hi <- eval_idx[eval_idx >= 2L * G]
    if (length(hi) == 0L) next
    D <- U0_full[hi + 1L, , drop = FALSE] - 2 * U0_full[hi - G + 1L, , drop = FALSE] + U0_full[hi - 2L * G + 1L, , drop = FALSE]
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

simulate_finite_alt_chunk <- function(seed_chunk, q_grid, m, k_max, gamma_vec, mosum_h_vec, weighted_cusum_omega_names,
                                      multiscale_h_sets, multiscale_scale_names, page_length_grid_size,
                                      weighted_length_grid_size, finite_eval_grid_size, ridge, range_floor,
                                      exact_page_scan, exact_weighted_scan) {
  lapply(seed_chunk, function(sd) {
    set.seed(sd)
    simulate_one_finite_alt_rep(q_grid = q_grid, m = m, k_max = k_max, gamma_vec = gamma_vec,
      mosum_h_vec = mosum_h_vec, weighted_cusum_omega_names = weighted_cusum_omega_names,
      multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names,
      page_length_grid_size = page_length_grid_size, weighted_length_grid_size = weighted_length_grid_size,
      finite_eval_grid_size = finite_eval_grid_size, ridge = ridge, range_floor = range_floor,
      exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan)
  })
}

simulate_finite_alt_critical_values_T <- function(T, q_grid, gamma_vec = c(0, 0.15), mosum_h_vec = c(0.10, 0.20),
                                                  weighted_cusum_omega_names = c("InvSqrt"),
                                                  multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                                  multiscale_scale_names = c("Equal"), alpha_levels = c(0.10, 0.05, 0.01),
                                                  nrep = 10000L, grid_monitor = 10000L, page_length_grid_size = 40L,
                                                  weighted_length_grid_size = 40L, finite_eval_grid_size = 250L,
                                                  ridge = 1e-10, range_floor = 1e-8, ncores = 1L,
                                                  seed = 123456789L, exact_page_scan = FALSE, exact_weighted_scan = FALSE,
                                                  batch_size = 100L, verbose = TRUE) {
  if (!is.finite(T)) stop("simulate_finite_alt_critical_values_T() requires finite T.")
  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  mosum_h_vec <- sort(unique(as.numeric(mosum_h_vec)))
  weighted_cusum_omega_names <- vapply(weighted_cusum_omega_names, normalize_omega_name, character(1L))
  multiscale_scale_names <- vapply(multiscale_scale_names, normalize_scale_weight_name, character(1L))
  alpha_levels <- as.numeric(alpha_levels)

  meta <- build_finite_alt_meta(q_grid, gamma_vec, mosum_h_vec, weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names)
  m <- compute_m_from_T(grid_monitor = grid_monitor, T = T)
  k_max <- as.integer(round(m * T))
  seeds <- make_seed_stream(nrep = nrep, seed = seed)
  core_plan <- resolve_safe_ncores(requested = ncores)
  ncores_use <- min(core_plan$used, batch_size)
  batches <- split(seeds, ceiling(seq_along(seeds) / batch_size))
  res_parts <- vector("list", length(batches))

  if (verbose) message(sprintf("Starting finite alternative-detector simulation: T=%s | q-grid=%s | nrep=%d | m=%d | k_max=%d | workers(requested=%d, used=%d) | batches=%d", canonical_T_scalar(T), format_q_grid(q_grid), nrep, m, k_max, core_plan$requested, ncores_use, length(batches)))

  for (b in seq_along(batches)) {
    seed_batch <- batches[[b]]
    if (ncores_use <= 1L || length(seed_batch) == 1L) {
      res_batch <- simulate_finite_alt_chunk(seed_batch, q_grid, m, k_max, gamma_vec, mosum_h_vec,
        weighted_cusum_omega_names, multiscale_h_sets, multiscale_scale_names, page_length_grid_size,
        weighted_length_grid_size, finite_eval_grid_size, ridge, range_floor, exact_page_scan, exact_weighted_scan)
    } else {
      cl <- make_psock_cluster(min(ncores_use, length(seed_batch)))
      on.exit(stop_psock_cluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c("simulate_one_finite_alt_rep", "simulate_finite_alt_chunk", "compute_metric_norms", "empty_value_store",
        "row_cumsum_weighted_squares", "leading_mahalanobis_sequence", "apply_scale_and_max", "col_maxs_base",
        "build_eval_index_grid", "build_linear_index_grid", "build_geometric_index_grid", "full_or_geometric_lag_grid", "g_gamma_sq",
        "make_weighted_cusum_omega", "make_multiscale_weight", "normalize_omega_name", "normalize_scale_weight_name",
        "q_grid", "m", "k_max", "gamma_vec", "mosum_h_vec", "weighted_cusum_omega_names", "multiscale_h_sets",
        "multiscale_scale_names", "page_length_grid_size", "weighted_length_grid_size", "finite_eval_grid_size", "ridge",
        "range_floor", "exact_page_scan", "exact_weighted_scan"), envir = environment())
      res_batch <- parallel::parLapplyLB(cl, seed_batch, function(sd) {
        set.seed(sd)
        simulate_one_finite_alt_rep(q_grid = q_grid, m = m, k_max = k_max, gamma_vec = gamma_vec,
          mosum_h_vec = mosum_h_vec, weighted_cusum_omega_names = weighted_cusum_omega_names,
          multiscale_h_sets = multiscale_h_sets, multiscale_scale_names = multiscale_scale_names,
          page_length_grid_size = page_length_grid_size, weighted_length_grid_size = weighted_length_grid_size,
          finite_eval_grid_size = finite_eval_grid_size, ridge = ridge, range_floor = range_floor,
          exact_page_scan = exact_page_scan, exact_weighted_scan = exact_weighted_scan)
      })
      stop_psock_cluster(cl)
      on.exit(NULL, add = FALSE)
    }
    res_parts[[b]] <- do.call(rbind, res_batch)
    if (verbose) message(sprintf("  T=%s batch %d/%d finished (%d/%d reps)", canonical_T_scalar(T), b, length(batches), min(b * batch_size, nrep), nrep))
  }

  res_mat <- do.call(rbind, res_parts)
  source_tag <- sprintf("simulated_alt_detectors_finite_T_%s_q_%d_%d_nrep_%d.csv", canonical_T_scalar(T), min(q_grid), max(q_grid), nrep)
  rows <- list(); idx <- 1L
  for (a in alpha_levels) {
    for (j in seq_len(nrow(meta))) {
      rowj <- meta[j, , drop = FALSE]
      rows[[idx]] <- data.frame(standardizer = rowj$standardizer, detector = rowj$detector, T = canonical_T_scalar(T), gamma = rowj$gamma,
        q = rowj$q, alpha = a, critical_value = quantile_upper(res_mat[, j], a), source_file = source_tag,
        bandwidth_h = rowj$bandwidth_h, omega_name = rowj$omega_name, hset_name = rowj$hset_name,
        scale_weight_name = rowj$scale_weight_name, exact_scan = FALSE, stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  if (verbose) message(sprintf("Finished finite alternative-detector simulation for T=%s.", canonical_T_scalar(T)))
  list(meta = meta, rows = if (length(rows)) do.call(rbind, rows) else empty_alt_table(), info = list(T = T, q_grid = q_grid, nrep = nrep, m = m, k_max = k_max, workers_used = ncores_use))
}
