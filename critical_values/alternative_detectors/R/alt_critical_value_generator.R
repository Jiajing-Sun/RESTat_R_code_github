# ==============================================================
# alt_critical_value_generator.R -- regenerate benchmark-detector
# critical values from Brownian-limit simulations
# ==============================================================

build_alt_cv_designs <- function(T_grid,
                                 q_grid,
                                 finite_nrep = 10000L,
                                 open_nrep = 5000L,
                                 finite_kmax = 10000L,
                                 open_train_grid_size = 1500L,
                                 open_monitor_grid_size = 2000L,
                                 open_s_max = 20) {
  q_grid <- sort(unique(as.integer(q_grid)))
  designs <- vector("list", length(T_grid))
  for (i in seq_along(T_grid)) {
    T_value <- T_grid[i]
    T_label <- canonical_T_scalar(T_value)
    if (identical(T_label, "Inf")) {
      r_grid <- seq_len(open_train_grid_size) / open_train_grid_size
      s_grid <- seq(open_s_max / open_monitor_grid_size, open_s_max, length.out = open_monitor_grid_size)
      designs[[i]] <- list(
        T_value = Inf,
        T_label = T_label,
        q_grid = q_grid,
        nrep = as.integer(open_nrep),
        r_grid = r_grid,
        s_grid = s_grid,
        source_file = sprintf("simulated_openend_q_%d_%d_nrep_%d.csv", min(q_grid), max(q_grid), as.integer(open_nrep))
      )
    } else {
      T_num <- as.numeric(T_value)
      m_grid_size <- as.integer(round(finite_kmax / T_num))
      if (abs(finite_kmax / T_num - m_grid_size) > 1e-10) {
        stop("finite_kmax / T must be integer for finite-horizon benchmark critical values.")
      }
      r_grid <- seq_len(m_grid_size) / m_grid_size
      s_grid <- seq_len(finite_kmax) / m_grid_size
      designs[[i]] <- list(
        T_value = T_num,
        T_label = T_label,
        q_grid = q_grid,
        nrep = as.integer(finite_nrep),
        r_grid = r_grid,
        s_grid = s_grid,
        source_file = sprintf("simulated_alt_detectors_finite_T_%s_q_%d_%d_nrep_%d.csv",
                              T_label, min(q_grid), max(q_grid), as.integer(finite_nrep))
      )
    }
  }
  designs
}

simulate_brownian_limit_bundle <- function(r_grid, s_grid, q_max) {
  times <- c(r_grid, 1 + s_grid)
  dtime <- diff(c(0, times))
  increments <- matrix(stats::rnorm(length(times) * q_max), nrow = length(times), ncol = q_max)
  increments <- sweep(increments, 1L, sqrt(dtime), FUN = "*")
  brownian <- apply(increments, 2L, cumsum)
  brownian <- safe_matrix(brownian, ncol = q_max)

  n_train <- length(r_grid)
  brownian_train <- brownian[seq_len(n_train), , drop = FALSE]
  brownian_monitor <- brownian[(n_train + 1L):nrow(brownian), , drop = FALSE]
  brownian_one <- brownian_train[n_train, ]
  bridge_train <- brownian_train - outer(r_grid, brownian_one)
  monitor_process <- brownian_monitor - outer(1 + s_grid, brownian_one)

  list(
    bridge_train = bridge_train,
    monitor_process = monitor_process,
    monitor_with_zero = rbind(rep(0, q_max), monitor_process),
    s_grid = s_grid
  )
}

build_alt_cv_contexts <- function(path_bundle, ridge = 1e-10, range_floor = 1e-8) {
  q_max <- ncol(path_bundle$monitor_process)
  bridge_train <- path_bundle$bridge_train
  d_ssms <- crossprod(bridge_train) / nrow(bridge_train)
  range_diag <- apply(bridge_train, 2L, max) - apply(bridge_train, 2L, min)
  range_diag[range_diag < range_floor] <- range_floor

  make_ctx <- function(standardizer, metric) {
    list(
      standardizer = standardizer,
      metric = metric,
      process = list(
        s = path_bundle$s_grid,
        k_max = length(path_bundle$s_grid),
        U0_full = path_bundle$monitor_with_zero
      )
    )
  }

  list(
    HAC = make_ctx("HAC", list(kind = "matrix", matrix = diag(q_max), ridge = ridge)),
    SSMS = make_ctx("SSMS", list(kind = "matrix", matrix = d_ssms, ridge = ridge)),
    RSMS = make_ctx("RSMS", list(kind = "diag", inv_diag = 1 / (range_diag^2)))
  )
}

quad_metric_leading <- function(D, context) {
  D <- safe_matrix(D)
  if (context$metric$kind == "matrix") {
    return(leading_mahalanobis_sequence(D, context$metric$matrix, ridge = context$metric$ridge))
  }
  row_cumsum_weighted_squares(D, w = context$metric$inv_diag)
}

grid_step_size <- function(s_grid) {
  if (length(s_grid) < 2L) stop("s_grid must contain at least two points.")
  stats::median(diff(s_grid))
}

bandwidth_to_steps <- function(h, s_grid) {
  ds <- grid_step_size(s_grid)
  g <- as.integer(round(h / ds))
  if (g < 1L) stop("Bandwidth ", h, " is too small for the monitoring grid.")
  if (abs(g * ds - h) > max(1e-8, ds * 1e-6)) {
    stop("Bandwidth ", h, " is not aligned with the monitoring grid step ", ds)
  }
  g
}

compute_alt_stat_leading <- function(method_row,
                                     context,
                                     page_length_grid_size = 40L,
                                     weighted_length_grid_size = 40L,
                                     exact_page_scan = FALSE,
                                     exact_weighted_scan = FALSE,
                                     multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20))) {
  proc <- context$process
  s <- proc$s
  s_full <- c(0, s)
  U0 <- proc$U0_full
  q_max <- ncol(U0)
  eval_idx <- seq_len(proc$k_max)
  detector <- method_row$detector

  if (detector == "FullCUSUM") {
    qmat <- quad_metric_leading(U0[eval_idx + 1L, , drop = FALSE], context)
    cur <- sweep(qmat, 1L, g_gamma_sq(s[eval_idx], method_row$gamma), FUN = "/")
    return(apply(cur, 2L, max))
  }

  if (detector == "PageCUSUM") {
    lag_grid <- full_or_geometric_lag_grid(proc$k_max, page_length_grid_size, exact = exact_page_scan)
    best <- rep(-Inf, q_max)
    for (L in lag_grid) {
      hi <- eval_idx[eval_idx >= L]
      if (length(hi) == 0L) next
      lo <- hi - L
      D <- U0[hi + 1L, , drop = FALSE] - U0[lo + 1L, , drop = FALSE]
      lens <- s_full[hi + 1L] - s_full[lo + 1L]
      qmat <- quad_metric_leading(D, context)
      cur <- sweep(qmat, 1L, g_gamma_sq(lens, method_row$gamma), FUN = "/")
      best <- pmax(best, apply(cur, 2L, max))
    }
    return(best)
  }

  if (detector == "WeightedCUSUM") {
    lag_grid <- full_or_geometric_lag_grid(proc$k_max, weighted_length_grid_size, exact = exact_weighted_scan)
    best <- rep(-Inf, q_max)
    for (L in lag_grid) {
      hi <- eval_idx[eval_idx >= L]
      if (length(hi) == 0L) next
      lo <- hi - L
      D <- U0[hi + 1L, , drop = FALSE] - U0[lo + 1L, , drop = FALSE]
      lens <- s_full[hi + 1L] - s_full[lo + 1L]
      qmat <- quad_metric_leading(D, context)
      cur <- sweep(qmat, 1L, g_gamma_sq(lens, method_row$gamma), FUN = "/")
      cur <- sweep(cur, 1L, make_weighted_cusum_omega(lens, method_row$omega_name), FUN = "*")
      best <- pmax(best, apply(cur, 2L, max))
    }
    return(best)
  }

  if (detector == "MOSUM") {
    h <- as.numeric(method_row$bandwidth_h)
    G <- bandwidth_to_steps(h, s)
    valid_hi <- eval_idx[eval_idx >= 2L * G]
    D <- U0[valid_hi + 1L, , drop = FALSE] -
      2 * U0[valid_hi - G + 1L, , drop = FALSE] +
      U0[valid_hi - 2L * G + 1L, , drop = FALSE]
    qmat <- quad_metric_leading(D, context)
    cur <- qmat / h
    return(apply(cur, 2L, max))
  }

  if (detector == "MultiscaleMOSUM") {
    h_vec <- multiscale_h_sets[[method_row$hset_name]]
    if (is.null(h_vec)) stop("Unknown multiscale h-set: ", method_row$hset_name)
    best <- rep(-Inf, q_max)
    for (h in h_vec) {
      G <- bandwidth_to_steps(h, s)
      valid_hi <- eval_idx[eval_idx >= 2L * G]
      D <- U0[valid_hi + 1L, , drop = FALSE] -
        2 * U0[valid_hi - G + 1L, , drop = FALSE] +
        U0[valid_hi - 2L * G + 1L, , drop = FALSE]
      qmat <- quad_metric_leading(D, context)
      cur <- qmat / h
      cur <- cur * make_multiscale_weight(h, method_row$scale_weight_name)
      best <- pmax(best, apply(cur, 2L, max))
    }
    return(best)
  }

  stop("Unknown benchmark detector: ", detector)
}

build_alt_stat_column_spec <- function(alt_catalog, q_grid) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(alt_catalog))) {
    row <- alt_catalog[i, , drop = FALSE]
    for (q in q_grid) {
      rows[[idx]] <- data.frame(
        standardizer = row$standardizer,
        detector = row$detector,
        gamma = row$gamma,
        q = as.integer(q),
        bandwidth_h = row$bandwidth_h,
        omega_name = row$omega_name,
        hset_name = row$hset_name,
        scale_weight_name = row$scale_weight_name,
        method_group = row$method_group,
        method_label = row$method_label,
        method_id = row$method_id,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$column_name <- sprintf("%s__q%02d", out$method_id, out$q)
  out
}

simulate_alt_stat_chunk <- function(seed_chunk,
                                    design,
                                    alt_catalog,
                                    q_grid,
                                    column_spec = NULL,
                                    ridge = 1e-10,
                                    range_floor = 1e-8,
                                    page_length_grid_size = 40L,
                                    weighted_length_grid_size = 40L,
                                    exact_page_scan = FALSE,
                                    exact_weighted_scan = FALSE,
                                    multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20))) {
  q_grid <- sort(unique(as.integer(q_grid)))
  q_max <- max(q_grid)
  if (is.null(column_spec)) column_spec <- build_alt_stat_column_spec(alt_catalog, q_grid)

  out <- matrix(NA_real_, nrow = length(seed_chunk), ncol = nrow(column_spec))
  colnames(out) <- column_spec$column_name

  for (rep_idx in seq_along(seed_chunk)) {
    set.seed(seed_chunk[rep_idx])
    path_bundle <- simulate_brownian_limit_bundle(design$r_grid, design$s_grid, q_max = q_max)
    contexts <- build_alt_cv_contexts(path_bundle, ridge = ridge, range_floor = range_floor)
    col_idx <- 1L
    for (i in seq_len(nrow(alt_catalog))) {
      row <- alt_catalog[i, , drop = FALSE]
      stats_q <- compute_alt_stat_leading(
        method_row = row,
        context = contexts[[row$standardizer]],
        page_length_grid_size = page_length_grid_size,
        weighted_length_grid_size = weighted_length_grid_size,
        exact_page_scan = exact_page_scan,
        exact_weighted_scan = exact_weighted_scan,
        multiscale_h_sets = multiscale_h_sets
      )
      out[rep_idx, col_idx:(col_idx + length(q_grid) - 1L)] <- stats_q[q_grid]
      col_idx <- col_idx + length(q_grid)
    }
  }

  out
}

simulate_alt_stat_matrix <- function(design,
                                     alt_catalog,
                                     q_grid,
                                     ridge = 1e-10,
                                     range_floor = 1e-8,
                                     page_length_grid_size = 40L,
                                     weighted_length_grid_size = 40L,
                                     exact_page_scan = FALSE,
                                     exact_weighted_scan = FALSE,
                                     multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                     progress_every = 100L,
                                     ncores = 1L,
                                     chunk_size = 100L,
                                     root = NULL) {
  q_grid <- sort(unique(as.integer(q_grid)))
  column_spec <- build_alt_stat_column_spec(alt_catalog, q_grid)
  chunk_size <- max(1L, as.integer(chunk_size))
  ncores <- max(1L, as.integer(ncores))
  seed_base <- as.integer(20260401L + sum(utf8ToInt(design$T_label)) * 100L)
  seeds <- make_seed_stream(design$nrep, seed = seed_base)
  seed_chunks <- split(seeds, ceiling(seq_along(seeds) / chunk_size))

  message(sprintf("Simulating benchmark critical values for T=%s with %d replications.", design$T_label, design$nrep))
  if (ncores <= 1L || length(seed_chunks) == 1L) {
    chunk_results <- vector("list", length(seed_chunks))
    completed <- 0L
    for (i in seq_along(seed_chunks)) {
      chunk_results[[i]] <- simulate_alt_stat_chunk(
        seed_chunk = seed_chunks[[i]],
        design = design,
        alt_catalog = alt_catalog,
        q_grid = q_grid,
        column_spec = column_spec,
        ridge = ridge,
        range_floor = range_floor,
        page_length_grid_size = page_length_grid_size,
        weighted_length_grid_size = weighted_length_grid_size,
        exact_page_scan = exact_page_scan,
        exact_weighted_scan = exact_weighted_scan,
        multiscale_h_sets = multiscale_h_sets
      )
      completed <- completed + length(seed_chunks[[i]])
      if (progress_every > 0L && (completed %% progress_every) == 0L) {
        message(sprintf("  completed %d / %d replications for T=%s", completed, design$nrep, design$T_label))
      }
    }
  } else {
    if (is.null(root) || !dir.exists(root)) {
      stop("A valid project root is required for parallel benchmark critical-value generation.")
    }

    cluster_error <- NULL
    cl <- tryCatch(
      make_psock_cluster(min(ncores, length(seed_chunks))),
      error = function(e) {
        cluster_error <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(cl)) {
      message("PSOCK cluster unavailable while regenerating benchmark critical values; falling back to serial execution. Details: ",
              cluster_error)
      chunk_results <- vector("list", length(seed_chunks))
      completed <- 0L
      for (i in seq_along(seed_chunks)) {
        chunk_results[[i]] <- simulate_alt_stat_chunk(
          seed_chunk = seed_chunks[[i]],
          design = design,
          alt_catalog = alt_catalog,
          q_grid = q_grid,
          column_spec = column_spec,
          ridge = ridge,
          range_floor = range_floor,
          page_length_grid_size = page_length_grid_size,
          weighted_length_grid_size = weighted_length_grid_size,
          exact_page_scan = exact_page_scan,
          exact_weighted_scan = exact_weighted_scan,
          multiscale_h_sets = multiscale_h_sets
        )
        completed <- completed + length(seed_chunks[[i]])
        if (progress_every > 0L && (completed %% progress_every) == 0L) {
          message(sprintf("  completed %d / %d replications for T=%s", completed, design$nrep, design$T_label))
        }
      }
    } else {
      on.exit(stop_psock_cluster(cl), add = TRUE)
      ROOT2 <- root
      parallel::clusterExport(cl, varlist = c("ROOT2"), envir = environment())
      parallel::clusterEvalQ(cl, {
        source(file.path(ROOT2, "R", "project_paths.R"), local = FALSE)
        source(file.path(ROOT2, "R", "utils.R"), local = FALSE)
        source(file.path(ROOT2, "R", "method_catalog.R"), local = FALSE)
        source(file.path(ROOT2, "R", "alt_critical_value_generator.R"), local = FALSE)
        NULL
      })

      chunk_results <- parallel::parLapplyLB(
        cl,
        seed_chunks,
        fun = function(sd_chunk,
                       design,
                       alt_catalog,
                       q_grid,
                       column_spec,
                       ridge,
                       range_floor,
                       page_length_grid_size,
                       weighted_length_grid_size,
                       exact_page_scan,
                       exact_weighted_scan,
                       multiscale_h_sets) {
          simulate_alt_stat_chunk(
            seed_chunk = sd_chunk,
            design = design,
            alt_catalog = alt_catalog,
            q_grid = q_grid,
            column_spec = column_spec,
            ridge = ridge,
            range_floor = range_floor,
            page_length_grid_size = page_length_grid_size,
            weighted_length_grid_size = weighted_length_grid_size,
            exact_page_scan = exact_page_scan,
            exact_weighted_scan = exact_weighted_scan,
            multiscale_h_sets = multiscale_h_sets
          )
        },
        design = design,
        alt_catalog = alt_catalog,
        q_grid = q_grid,
        column_spec = column_spec,
        ridge = ridge,
        range_floor = range_floor,
        page_length_grid_size = page_length_grid_size,
        weighted_length_grid_size = weighted_length_grid_size,
        exact_page_scan = exact_page_scan,
        exact_weighted_scan = exact_weighted_scan,
        multiscale_h_sets = multiscale_h_sets
      )
      stop_psock_cluster(cl)
      on.exit(NULL, add = FALSE)
    }
  }

  out <- do.call(rbind, chunk_results)
  list(stats = out, column_spec = column_spec)
}

exact_scan_flag <- function(detector) {
  if (detector %in% c("PageCUSUM", "WeightedCUSUM")) return(FALSE)
  TRUE
}

summarize_alt_stat_matrix <- function(stat_matrix,
                                      column_spec,
                                      alpha_levels,
                                      T_label,
                                      source_file) {
  rows <- list()
  idx <- 1L
  probs <- 1 - sort(unique(alpha_levels), decreasing = FALSE)
  alpha_levels <- sort(unique(alpha_levels), decreasing = FALSE)

  for (j in seq_len(nrow(column_spec))) {
    qs <- as.numeric(stats::quantile(stat_matrix[, j], probs = probs, names = FALSE, type = 7, na.rm = TRUE))
    for (k in seq_along(alpha_levels)) {
      rows[[idx]] <- data.frame(
        standardizer = column_spec$standardizer[j],
        detector = column_spec$detector[j],
        T = T_label,
        gamma = column_spec$gamma[j],
        q = column_spec$q[j],
        alpha = alpha_levels[k],
        critical_value = qs[k],
        source_file = source_file,
        bandwidth_h = column_spec$bandwidth_h[j],
        omega_name = ifelse(column_spec$omega_name[j] == "", NA_character_, column_spec$omega_name[j]),
        hset_name = ifelse(column_spec$hset_name[j] == "", NA_character_, column_spec$hset_name[j]),
        scale_weight_name = ifelse(column_spec$scale_weight_name[j] == "", NA_character_, column_spec$scale_weight_name[j]),
        exact_scan = exact_scan_flag(column_spec$detector[j]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, rows)
  out[order(out$standardizer, out$detector, out$T, out$q, out$alpha, out$gamma), , drop = FALSE]
}

generate_alt_detector_critical_values <- function(T_grid = c(1, 2, 5, 10),
                                                  q_grid = 1:30,
                                                  alpha_levels = c(0.01, 0.05, 0.10),
                                                  gamma_vec = c(0, 0.15),
                                                  mosum_h_vec = c(0.10, 0.20),
                                                  weighted_omega_names = c("InvSqrt"),
                                                  multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                                  multiscale_scale_names = c("Equal"),
                                                  finite_nrep = 10000L,
                                                  open_nrep = 5000L,
                                                  finite_kmax = 10000L,
                                                  open_train_grid_size = 1500L,
                                                  open_monitor_grid_size = 2000L,
                                                  open_s_max = 20,
                                                  page_length_grid_size = 40L,
                                                  weighted_length_grid_size = 40L,
                                                  exact_page_scan = FALSE,
                                                  exact_weighted_scan = FALSE,
                                                  ridge = 1e-10,
                                                  range_floor = 1e-8,
                                                  progress_every = 100L,
                                                  ncores = 1L,
                                                  chunk_size = 100L,
                                                  root = NULL) {
  alt_catalog <- build_alt_method_catalog(
    gamma_vec = gamma_vec,
    mosum_h_vec = mosum_h_vec,
    weighted_omega_names = weighted_omega_names,
    multiscale_h_sets = multiscale_h_sets,
    multiscale_scale_names = multiscale_scale_names
  )
  designs <- build_alt_cv_designs(
    T_grid = T_grid,
    q_grid = q_grid,
    finite_nrep = finite_nrep,
    open_nrep = open_nrep,
    finite_kmax = finite_kmax,
    open_train_grid_size = open_train_grid_size,
    open_monitor_grid_size = open_monitor_grid_size,
    open_s_max = open_s_max
  )

  pieces <- lapply(designs, function(design) {
    sim <- simulate_alt_stat_matrix(
      design = design,
      alt_catalog = alt_catalog,
      q_grid = q_grid,
      ridge = ridge,
      range_floor = range_floor,
      page_length_grid_size = page_length_grid_size,
      weighted_length_grid_size = weighted_length_grid_size,
      exact_page_scan = exact_page_scan,
      exact_weighted_scan = exact_weighted_scan,
      multiscale_h_sets = multiscale_h_sets,
      progress_every = progress_every,
      ncores = ncores,
      chunk_size = chunk_size,
      root = root
    )
    summarize_alt_stat_matrix(
      stat_matrix = sim$stats,
      column_spec = sim$column_spec,
      alpha_levels = alpha_levels,
      T_label = design$T_label,
      source_file = design$source_file
    )
  })

  do.call(rbind, pieces)
}
