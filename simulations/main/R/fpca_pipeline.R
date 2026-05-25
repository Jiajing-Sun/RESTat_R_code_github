# ==============================================================
# fpca_pipeline.R -- training-sample FPCA and score construction
# ============================================================== 

fpca_project_scores <- function(fd_data,
                                m,
                                q_cap = 30L,
                                fve_threshold = 0.95,
                                fixed_q = NA_integer_,
                                n_grid = 301L) {
  t_grid <- seq(0, 1, length.out = n_grid)
  values_all <- t(fda::eval.fd(t_grid, fd_data))
  values_all <- safe_matrix(values_all)
  total_n <- nrow(values_all)
  if (m >= total_n) stop("m must be smaller than the total sample size.")

  values_train <- values_all[1:m, , drop = FALSE]
  values_monitor <- values_all[(m + 1L):total_n, , drop = FALSE]
  Ly_train <- split(values_train, row(values_train))
  Lt_train <- replicate(length(Ly_train), t_grid, simplify = FALSE)

  fpca_train <- fdapace::FPCA(
    Ly = Ly_train,
    Lt = Lt_train,
    optns = list(methodMuCovEst = "smooth", FVEthreshold = fve_threshold, methodSelectK = "FVE")
  )

  mean_train <- as.numeric(fpca_train$mu)
  phi_train <- safe_matrix(fpca_train$phi)
  scores_train <- safe_matrix(fpca_train$xiEst)
  centered_monitor <- sweep(values_monitor, 2, mean_train, FUN = "-")
  scores_monitor <- safe_matrix(centered_monitor %*% phi_train * (1 / length(t_grid)))

  q_available <- min(ncol(scores_train), ncol(phi_train), q_cap)
  if (is.na(fixed_q)) {
    q_used <- q_available
  } else {
    q_used <- min(as.integer(fixed_q[1L]), q_available)
  }
  if (q_used < 1L) stop("FPCA produced zero retained components.")

  scores_train <- scores_train[, 1:q_used, drop = FALSE]
  scores_monitor <- scores_monitor[, 1:q_used, drop = FALSE]

  list(
    t_grid = t_grid,
    values_all = values_all,
    values_train = values_train,
    values_monitor = values_monitor,
    mean_train = mean_train,
    phi_train = phi_train[, 1:q_used, drop = FALSE],
    scores_train = scores_train,
    scores_monitor = scores_monitor,
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
