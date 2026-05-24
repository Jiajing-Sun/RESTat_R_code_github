# ==============================================================
# openend_critical_values.R -- open-end (T = Inf) critical values
# ============================================================== 

build_openend_meta <- function(q_grid,
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

simulate_one_openend_rep <- function(q_grid,
                                     n_train_grid = 1500L,
                                     n_open_grid = 2000L,
                                     gamma_vec = c(0, 0.15),
                                     weight_names = c("U", "Early", "Mid", "Late"),
                                     ridge = 1e-10,
                                     range_floor = 1e-8) {
  q_grid <- sort(unique(as.integer(q_grid)))
  q_max <- max(q_grid)

  dB <- matrix(rnorm(n_train_grid * q_max, sd = 1 / sqrt(n_train_grid)),
               nrow = n_train_grid, ncol = q_max)
  B <- apply(dB, 2, cumsum)
  r <- (1:n_train_grid) / n_train_grid
  Z <- B[n_train_grid, ]
  B0 <- B - outer(r, Z)

  V <- crossprod(B0) / n_train_grid
  ranges <- apply(B0, 2, max) - apply(B0, 2, min)
  ranges[ranges < range_floor] <- range_floor
  inv2 <- 1 / (ranges^2)

  dG <- matrix(rnorm(n_open_grid * q_max, sd = 1 / sqrt(n_open_grid)),
               nrow = n_open_grid, ncol = q_max)
  G <- apply(dG, 2, cumsum)

  x <- c(0, (1:n_open_grid) / n_open_grid)
  G <- rbind(rep(0, q_max), G)

  w_mat <- sapply(weight_names, function(w) make_cvm_weight_open_x(x, weight = w))
  w_mat <- as.matrix(w_mat)
  colnames(w_mat) <- vapply(weight_names, normalize_weight_name, character(1L))

  values <- numeric(0L)

  for (q in q_grid) {
    Gq <- G[, 1:q, drop = FALSE]

    quad_ss <- rowSums((Gq %*% safe_solve(V[1:q, 1:q, drop = FALSE], ridge = ridge)) * Gq)
    quad_rs <- rowSums(sweep(Gq^2, 2, inv2[1:q], FUN = "*"))
    quad_hac <- rowSums(Gq^2)

    for (g in gamma_vec) {
      denom <- x^(2 * g)
      denom[1L] <- Inf
      values <- c(values,
                  max(quad_ss / denom),
                  max(quad_rs / denom),
                  max(quad_hac / denom))
    }

    for (j in seq_len(ncol(w_mat))) {
      wj <- w_mat[, j]
      values <- c(values,
                  trapz_equal_grid(wj * quad_ss),
                  trapz_equal_grid(wj * quad_rs),
                  trapz_equal_grid(wj * quad_hac))
    }
  }

  values
}

simulate_openend_chunk <- function(seed_chunk,
                                   q_grid,
                                   n_train_grid,
                                   n_open_grid,
                                   gamma_vec,
                                   weight_names,
                                   ridge,
                                   range_floor) {
  lapply(seed_chunk, function(sd) {
    set.seed(sd)
    simulate_one_openend_rep(
      q_grid = q_grid,
      n_train_grid = n_train_grid,
      n_open_grid = n_open_grid,
      gamma_vec = gamma_vec,
      weight_names = weight_names,
      ridge = ridge,
      range_floor = range_floor
    )
  })
}

simulate_openend_critical_values <- function(q_grid,
                                             gamma_vec = c(0, 0.15),
                                             weight_names = c("U", "Early", "Mid", "Late"),
                                             alpha_levels = c(0.10, 0.05, 0.01),
                                             nrep = 5000L,
                                             n_train_grid = 1500L,
                                             n_open_grid = 2000L,
                                             ridge = 1e-10,
                                             range_floor = 1e-8,
                                             ncores = 1L,
                                             seed = 13579L,
                                             verbose = TRUE) {
  q_grid <- sort(unique(as.integer(q_grid)))
  gamma_vec <- as.numeric(gamma_vec)
  weight_names <- vapply(weight_names, normalize_weight_name, character(1L))
  alpha_levels <- as.numeric(alpha_levels)
  meta <- build_openend_meta(q_grid = q_grid, gamma_vec = gamma_vec, weight_names = weight_names)
  seeds <- make_seed_stream(nrep = nrep, seed = seed)

  core_plan <- resolve_safe_ncores(requested = ncores)
  ncores_use <- min(core_plan$used, length(seeds))
  nchunks <- min(length(seeds), max(1L, 4L * ncores_use))
  seed_chunks <- split_into_chunks(seeds, nchunks = nchunks)

  if (verbose) {
    message(sprintf(
      paste0(
        "Starting open-end simulation: T=Inf | q-grid=%s | nrep=%d | train-grid=%d | open-grid=%d | ",
        "workers(requested=%d, used=%d) | chunks=%d"
      ),
      format_q_grid(q_grid), nrep, n_train_grid, n_open_grid,
      core_plan$requested, ncores_use, length(seed_chunks)
    ))
  }

  if (ncores_use <= 1L) {
    res_list <- simulate_openend_chunk(
      seed_chunk = seeds,
      q_grid = q_grid,
      n_train_grid = n_train_grid,
      n_open_grid = n_open_grid,
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
        "simulate_one_openend_rep", "simulate_openend_chunk", "safe_solve", "trapz_equal_grid", "make_cvm_weight_open_x", "normalize_weight_name",
        "q_grid", "n_train_grid", "n_open_grid", "gamma_vec", "weight_names", "ridge", "range_floor"
      ),
      envir = environment()
    )

    res_chunks <- parallel::parLapplyLB(cl, seed_chunks, function(sd_chunk) {
      simulate_openend_chunk(
        seed_chunk = sd_chunk,
        q_grid = q_grid,
        n_train_grid = n_train_grid,
        n_open_grid = n_open_grid,
        gamma_vec = gamma_vec,
        weight_names = weight_names,
        ridge = ridge,
        range_floor = range_floor
      )
    })
    res_list <- unlist(res_chunks, recursive = FALSE, use.names = FALSE)
  }

  res_mat <- do.call(rbind, res_list)
  source_tag <- sprintf("simulated_openend_q_%d_%d_nrep_%d.csv", min(q_grid), max(q_grid), nrep)

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
          T = "Inf",
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
          T = "Inf",
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
    message("Finished open-end simulation for T=Inf.")
  }

  list(
    meta = meta,
    base_rows = if (length(base_rows)) do.call(rbind, base_rows) else empty_base_table(),
    weight_rows = if (length(weight_rows)) do.call(rbind, weight_rows) else empty_weight_table(),
    info = list(T = Inf, q_grid = q_grid, nrep = nrep,
                n_train_grid = n_train_grid, n_open_grid = n_open_grid,
                workers_used = ncores_use)
  )
}
