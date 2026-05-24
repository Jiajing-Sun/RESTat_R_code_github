# ==============================================================
# fpca_pipeline.R -- empirical FPCA from discretized curves using
# base-R SVD on the frozen training sample
# ============================================================== 

fpca_project_scores <- function(curve_values,
                                m,
                                q_cap = 30L,
                                fve_threshold = 0.95,
                                fixed_q = NA_integer_,
                                n_grid = NULL) {
  values_all <- safe_matrix(curve_values)
  total_n <- nrow(values_all)
  if (m >= total_n) stop("m must be smaller than the total sample size.")
  n_grid_eff <- if (is.null(n_grid)) ncol(values_all) else min(as.integer(n_grid), ncol(values_all))
  values_all <- values_all[, seq_len(n_grid_eff), drop = FALSE]

  values_train <- values_all[1:m, , drop = FALSE]
  values_monitor <- values_all[(m + 1L):total_n, , drop = FALSE]

  mean_train <- colMeans(values_train)
  centered_train <- sweep(values_train, 2, mean_train, FUN = "-")
  centered_monitor <- sweep(values_monitor, 2, mean_train, FUN = "-")
  grid_scale <- sqrt(ncol(values_all))

  sv <- svd(centered_train / grid_scale, nu = 0L, nv = min(nrow(centered_train), ncol(centered_train), q_cap))
  phi_all <- safe_matrix(sv$v)
  scores_train_all <- safe_matrix((centered_train / grid_scale) %*% phi_all)
  scores_monitor_all <- safe_matrix((centered_monitor / grid_scale) %*% phi_all)
  var_all <- sv$d^2

  q_available <- min(ncol(scores_train_all), q_cap)
  if (q_available < 1L) stop("Empirical FPCA produced zero retained components.")

  if (is.na(fixed_q)) {
    cprop <- cumsum(var_all[seq_len(q_available)]) / sum(var_all[seq_len(q_available)])
    q_fve <- which(cprop >= fve_threshold)[1L]
    if (is.na(q_fve) || q_fve < 1L) q_fve <- q_available
    q_used <- min(q_fve, q_available)
  } else {
    q_used <- min(as.integer(fixed_q[1L]), q_available)
  }
  q_used <- max(1L, q_used)

  list(
    values_all = values_all,
    values_train = values_train,
    values_monitor = values_monitor,
    mean_train = mean_train,
    phi_train = phi_all[, 1:q_used, drop = FALSE],
    scores_train = scores_train_all[, 1:q_used, drop = FALSE],
    scores_monitor = scores_monitor_all[, 1:q_used, drop = FALSE],
    q_used = q_used,
    q_available = q_available
  )
}

whiten_scores <- function(scores_train, scores_monitor, ridge = 1e-8) {
  Zt <- safe_matrix(scores_train)
  Zm <- safe_matrix(scores_monitor)
  S <- stats::cov(Zt)
  S <- safe_matrix(S, ncol = ncol(Zt))
  S <- (S + t(S)) / 2 + ridge * diag(ncol(S))
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, ridge)
  inv_sqrt <- ee$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
  list(
    scores_train = Zt %*% inv_sqrt,
    scores_monitor = Zm %*% inv_sqrt,
    transform = inv_sqrt
  )
}
