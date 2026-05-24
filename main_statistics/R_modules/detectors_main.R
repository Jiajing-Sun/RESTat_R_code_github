# ==============================================================
# detectors_main.R -- main-paper KS and weighted-CvM statistics
# ============================================================== 

center_train_monitor <- function(scores_train, scores_monitor) {
  Zt <- safe_matrix(scores_train)
  Zm <- safe_matrix(scores_monitor)
  mu <- colMeans(Zt)
  list(
    train = sweep(Zt, 2, mu, FUN = "-"),
    monitor = sweep(Zm, 2, mu, FUN = "-")
  )
}

build_score_processes <- function(scores_train, scores_monitor) {
  Zt <- safe_matrix(scores_train)
  Zm <- safe_matrix(scores_monitor)
  m <- nrow(Zt)
  k_max <- nrow(Zm)
  q <- ncol(Zt)
  c_train <- safe_matrix(apply(Zt, 2, cumsum), ncol = q)
  c_monitor <- safe_matrix(apply(Zm, 2, cumsum), ncol = q)
  train_total <- colSums(Zt)
  r <- (1:m) / m
  s <- (1:k_max) / m
  B0_train <- c_train - outer(r, train_total)
  U_full <- c_monitor - outer(s, train_total)
  U0_full <- rbind(rep(0, q), U_full)
  list(m = m, k_max = k_max, q = q, c_train = c_train, c_monitor = c_monitor,
       train_total = train_total, r = r, s = s, B0_train = B0_train,
       U_full = U_full, U0_full = U0_full)
}

prepare_standardizer_context <- function(scores_train, scores_monitor, standardizer,
                                         ridge = 1e-10, range_floor = 1e-8,
                                         hac_bandwidth = NULL) {
  centered <- center_train_monitor(scores_train, scores_monitor)
  proc <- build_score_processes(centered$train, centered$monitor)
  std_u <- toupper(standardizer)
  if (std_u == "HAC") {
    Sigma <- hac_bartlett(centered$train, bandwidth = hac_bandwidth, center = FALSE, ridge = ridge)
    metric <- list(kind = "matrix", inv = safe_solve(Sigma, ridge = ridge))
  } else if (std_u == "SSMS") {
    D <- crossprod(proc$B0_train) / (proc$m^2)
    metric <- list(kind = "matrix", inv = safe_solve(D, ridge = ridge))
  } else if (std_u == "RSMS") {
    rg <- apply(proc$B0_train, 2, max) - apply(proc$B0_train, 2, min)
    rg[rg < range_floor] <- range_floor
    metric <- list(kind = "diag", inv_diag = 1 / (rg^2))
  } else {
    stop("Unknown standardizer: ", standardizer)
  }
  list(standardizer = standardizer, process = proc, metric = metric)
}

quad_metric_rows <- function(D, context) {
  D <- safe_matrix(D)
  std_u <- toupper(context$standardizer)
  if (context$metric$kind == "matrix") {
    qv <- rowSums((D %*% context$metric$inv) * D) / context$process$m
  } else {
    qv <- rowSums(sweep(D^2, 2, context$metric$inv_diag, FUN = "*"))
  }
  as.numeric(qv)
}

evaluate_main_method <- function(method_row, context_map, q_used, T_value, alpha_levels, cv_main) {
  ctx <- context_map[[toupper(method_row$standardizer)]]
  s <- ctx$process$s
  k_index <- seq_len(ctx$process$k_max)
  quad <- quad_metric_rows(ctx$process$U_full[, 1:q_used, drop = FALSE], ctx)

  if (toupper(method_row$type) == "KS") {
    path <- quad / g_gamma_sq(s, method_row$gamma)
    stat <- max(path)
  } else {
    w <- make_cvm_weight_finite(k_index, m = ctx$process$m, T = T_value, weight = method_row$weight_name)
    path <- cumsum(w * (quad / ((1 + s)^2))) / ctx$process$m
    stat <- max(path)
  }

  out <- vector("list", length(alpha_levels))
  for (i in seq_along(alpha_levels)) {
    a <- alpha_levels[i]
    cv <- lookup_main_critical_value(cv_main, standardizer = method_row$standardizer, type = method_row$type,
                                     T = T_value, q = q_used, alpha = a, gamma = method_row$gamma,
                                     weight_name = method_row$weight_name)
    out[[i]] <- data.frame(
      family = method_row$family,
      standardizer = method_row$standardizer,
      detector = method_row$detector,
      type = method_row$type,
      gamma = method_row$gamma,
      weight_name = ifelse(method_row$weight_name == "", NA, method_row$weight_name),
      bandwidth_h = NA_real_,
      omega_name = NA_character_,
      hset_name = NA_character_,
      scale_weight_name = NA_character_,
      method_group = method_row$method_group,
      method_label = method_row$method_label,
      method_id = method_row$method_id,
      alpha = a,
      statistic = stat,
      critical_value = cv,
      reject = stat > cv,
      first_rejection = first_crossing(path, cv, k_index = k_index),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
