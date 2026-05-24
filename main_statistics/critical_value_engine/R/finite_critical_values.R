# ==============================================================
# finite_critical_values.R -- finite-horizon critical values
# ============================================================== 

compute_m_from_T <- function(grid_monitor = 10000L, T) {
  if (!is.finite(T)) stop("compute_m_from_T() is for finite T only.")
  m <- grid_monitor / T
  if (abs(m - round(m)) > 1e-8) {
    stop(sprintf("grid_monitor=%s is not divisible by T=%s.", grid_monitor, T))
  }
  as.integer(round(m))
}

build_finite_meta <- function(q_grid,
                              gamma_vec,
                              weight_names = c("U", "Early", "Mid", "Late")) {
  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  weight_names <- vapply(weight_names, normalize_weight_name, character(1L))

  rows <- vector("list", length = 3L * length(q_grid) * length(gamma_vec) + 3L * length(q_grid) * length(weight_names))
  idx <- 1L

  for (q in q_grid) {
    for (g in gamma_vec) {
      rows[[idx]] <- data.frame(stat = "SSMS", type = "KS",  gamma = g, q = q, weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(stat = "RSMS", type = "KS",  gamma = g, q = q, weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(stat = "HAC",  type = "KS",  gamma = g, q = q, weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
    for (w in weight_names) {
      rows[[idx]] <- data.frame(stat = "SSMS", type = "CvM", gamma = 0, q = q, weight_name = w, stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(stat = "RSMS", type = "CvM", gamma = 0, q = q, weight_name = w, stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(stat = "HAC",  type = "CvM", gamma = 0, q = q, weight_name = w, stringsAsFactors = FALSE); idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

simulate_one_finite_rep <- function(q_grid,
                                    m,
                                    k_max,
                                    gamma_vec,
                                    weight_names,
                                    ridge = 1e-10,
                                    range_floor = 1e-8) {
  q_grid <- sort(unique(as.integer(q_grid)))
  q_max <- max(q_grid)

  total_steps <- m + k_max
  dW <- matrix(rnorm(total_steps * q_max), nrow = total_steps, ncol = q_max) / sqrt(m)
  W  <- apply(dW, 2, cumsum)

  r <- (1:m) / m
  s <- (1:k_max) / m
  one_plus_s <- 1 + s
  ratio <- s / one_plus_s

  B1 <- W[m, ]
  B0 <- W[1:m, , drop = FALSE] - outer(r, B1)

  V <- crossprod(B0) / m
  ranges <- apply(B0, 2, max) - apply(B0, 2, min)
  ranges[ranges < range_floor] <- range_floor
  inv2 <- 1 / (ranges^2)

  W_mon <- W[(m + 1):(m + k_max), , drop = FALSE]
  U <- W_mon - outer(one_plus_s, B1)

  denom_ks <- lapply(gamma_vec, function(g) (one_plus_s)^2 * (ratio)^(2 * g))
  denom_cvm <- (one_plus_s)^2

  w_mat <- sapply(weight_names, function(w) make_cvm_weight_finite(1:k_max, m = m, T = k_max / m, weight = w))
  w_mat <- as.matrix(w_mat)
  colnames(w_mat) <- vapply(weight_names, normalize_weight_name, character(1L))

  values <- numeric(0L)

  for (q in q_grid) {
    Uq <- U[, 1:q, drop = FALSE]
    Vq_inv <- safe_solve(V[1:q, 1:q, drop = FALSE], ridge = ridge)

    quad_ss <- rowSums((Uq %*% Vq_inv) * Uq)
    quad_rs <- rowSums(sweep(Uq^2, 2, inv2[1:q], FUN = "*"))
    quad_hac <- rowSums(Uq^2)

    for (gi in seq_along(gamma_vec)) {
      dgi <- denom_ks[[gi]]
      values <- c(values,
                  max(quad_ss / dgi),
                  max(quad_rs / dgi),
                  max(quad_hac / dgi))
    }

    M_ss0 <- quad_ss / denom_cvm
    M_rs0 <- quad_rs / denom_cvm
    M_hac0 <- quad_hac / denom_cvm

    for (j in seq_len(ncol(w_mat))) {
      wj <- w_mat[, j]
      I_ss <- cumsum(wj * M_ss0) / m
      I_rs <- cumsum(wj * M_rs0) / m
      I_hac <- cumsum(wj * M_hac0) / m
      values <- c(values,
                  max(I_ss),
                  max(I_rs),
                  max(I_hac))
    }
  }

  values
}

simulate_finite_chunk <- function(seed_chunk,
                                  q_grid,
                                  m,
                                  k_max,
                                  gamma_vec,
                                  weight_names,
                                  ridge,
                                  range_floor) {
  lapply(seed_chunk, function(sd) {
    set.seed(sd)
    simulate_one_finite_rep(
      q_grid = q_grid,
      m = m,
      k_max = k_max,
      gamma_vec = gamma_vec,
      weight_names = weight_names,
      ridge = ridge,
      range_floor = range_floor
    )
  })
}

simulate_finite_critical_values_T <- function(T,
                                              q_grid,
                                              gamma_vec = c(0, 0.15),
                                              weight_names = c("U", "Early", "Mid", "Late"),
                                              alpha_levels = c(0.10, 0.05, 0.01),
                                              nrep = 10000L,
                                              grid_monitor = 10000L,
                                              ridge = 1e-10,
                                              range_floor = 1e-8,
                                              ncores = 1L,
                                              seed = 123456789L,
                                              verbose = TRUE) {
  if (!is.finite(T)) stop("simulate_finite_critical_values_T() requires finite T.")

  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  weight_names <- vapply(weight_names, normalize_weight_name, character(1L))
  alpha_levels <- as.numeric(alpha_levels)

  meta <- build_finite_meta(q_grid = q_grid, gamma_vec = gamma_vec, weight_names = weight_names)

  m <- compute_m_from_T(grid_monitor = grid_monitor, T = T)
  k_max <- as.integer(round(m * T))
  seeds <- make_seed_stream(nrep = nrep, seed = seed)

  core_plan <- resolve_safe_ncores(requested = ncores)
  ncores_use <- min(core_plan$used, length(seeds))
  nchunks <- min(length(seeds), max(1L, 4L * ncores_use))
  seed_chunks <- split_into_chunks(seeds, nchunks = nchunks)

  if (verbose) {
    message(sprintf(
      paste0(
        "Starting finite-horizon simulation: T=%s | q-grid=%s | nrep=%d | m=%d | k_max=%d | ",
        "workers(requested=%d, used=%d) | chunks=%d"
      ),
      canonical_T_scalar(T), format_q_grid(q_grid), nrep, m, k_max,
      core_plan$requested, ncores_use, length(seed_chunks)
    ))
  }

  if (ncores_use <= 1L) {
    res_list <- simulate_finite_chunk(
      seed_chunk = seeds,
      q_grid = q_grid,
      m = m,
      k_max = k_max,
      gamma_vec = gamma_vec,
      weight_names = weight_names,
      ridge = ridge,
      range_floor = range_floor
    )
  } else {
    cl <- make_psock_cluster(ncores_use)
    on.exit(stop_psock_cluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "simulate_one_finite_rep", "simulate_finite_chunk", "safe_solve", "make_cvm_weight_finite", "normalize_weight_name",
        "q_grid", "m", "k_max", "gamma_vec", "weight_names", "ridge", "range_floor"
      ),
      envir = environment()
    )

    res_chunks <- parallel::parLapplyLB(cl, seed_chunks, function(sd_chunk) {
      simulate_finite_chunk(
        seed_chunk = sd_chunk,
        q_grid = q_grid,
        m = m,
        k_max = k_max,
        gamma_vec = gamma_vec,
        weight_names = weight_names,
        ridge = ridge,
        range_floor = range_floor
      )
    })
    res_list <- unlist(res_chunks, recursive = FALSE, use.names = FALSE)
  }

  res_mat <- do.call(rbind, res_list)
  source_tag <- sprintf("simulated_finite_T_%s_q_%d_%d_nrep_%d.csv",
                        canonical_T_scalar(T), min(q_grid), max(q_grid), nrep)

  base_rows <- list()
  weight_rows <- list()
  ib <- 1L
  iw <- 1L

  for (a in alpha_levels) {
    for (j in seq_len(nrow(meta))) {
      rowj <- meta[j, , drop = FALSE]
      qv <- quantile_upper(res_mat[, j], a)

      if (toupper(rowj$type) == "KS" || rowj$weight_name == "U") {
        base_rows[[ib]] <- data.frame(
          stat = rowj$stat,
          type = rowj$type,
          T = canonical_T_scalar(T),
          gamma = rowj$gamma,
          q = rowj$q,
          alpha = a,
          critical_value = qv,
          source_file = source_tag,
          stringsAsFactors = FALSE
        )
        ib <- ib + 1L
      } else {
        weight_rows[[iw]] <- data.frame(
          stat = rowj$stat,
          type = rowj$type,
          T = canonical_T_scalar(T),
          gamma = rowj$gamma,
          q = rowj$q,
          alpha = a,
          critical_value = qv,
          source_file = source_tag,
          weight_name = rowj$weight_name,
          stringsAsFactors = FALSE
        )
        iw <- iw + 1L
      }
    }
  }

  if (verbose) {
    message(sprintf("Finished finite-horizon simulation for T=%s.", canonical_T_scalar(T)))
  }

  list(
    meta = meta,
    base_rows = if (length(base_rows)) do.call(rbind, base_rows) else empty_base_table(),
    weight_rows = if (length(weight_rows)) do.call(rbind, weight_rows) else empty_weight_table(),
    info = list(T = T, m = m, k_max = k_max, q_grid = q_grid, nrep = nrep, workers_used = ncores_use)
  )
}
