# ==============================================================
# genData.R -- matrix-based functional data generators
# No fda/fdapace dependency: this is deliberate for Windows PSOCK
# robustness.
# ============================================================== 

make_basis_matrix <- function(t_grid, nbasis = 21L, basis_type = c("bspline", "fourier")) {
  basis_type <- match.arg(basis_type)
  t_grid <- as.numeric(t_grid)
  nbasis <- as.integer(nbasis)
  if (basis_type == "bspline") {
    B <- splines::bs(t_grid, df = nbasis, degree = 3L, intercept = TRUE)
  } else {
    j <- seq_len(nbasis)
    B <- sapply(j, function(k) {
      if (k == 1L) rep(1, length(t_grid)) else sqrt(2) * sin(pi * (k - 1L) * t_grid)
    })
  }
  B <- as.matrix(B)
  qrB <- qr(B)
  Q <- qr.Q(qrB)
  Q[, seq_len(min(nbasis, ncol(Q))), drop = FALSE]
}

make_mu_values <- function(t_grid, muInfo) {
  mu_type <- muInfo$type %||% "sin"
  mu_a <- as.numeric(muInfo$a %||% 0)
  if (identical(mu_a, 0)) return(rep(0, length(t_grid)))
  if (mu_type == "sin") return(mu_a * sin(2 * pi * t_grid))
  if (mu_type == "horv") return(mu_a * t_grid * (1 - t_grid))
  if (mu_type == "const") return(rep(mu_a, length(t_grid)))
  stop("Unsupported muInfo$type: ", mu_type)
}

generate_bb_curves <- function(total_n, t_grid, mean_curve = NULL) {
  n_grid <- length(t_grid)
  dt <- diff(t_grid)
  eps <- matrix(stats::rnorm(total_n * (n_grid - 1L)), nrow = total_n, ncol = n_grid - 1L)
  eps <- sweep(eps, 2, sqrt(pmax(dt, 1e-8)), FUN = "*")
  W <- t(apply(eps, 1, cumsum))
  W <- cbind(0, W)
  bridge <- W - outer(W[, n_grid], t_grid)
  if (!is.null(mean_curve)) bridge <- sweep(bridge, 2, mean_curve, FUN = "+")
  bridge
}

generate_iid_curves <- function(total_n, basis_matrix, sd_pattern, mean_curve = NULL, gaussian = TRUE) {
  nbasis <- ncol(basis_matrix)
  if (gaussian) {
    coef <- matrix(stats::rnorm(total_n * nbasis), nrow = total_n, ncol = nbasis)
  } else {
    coef <- matrix(stats::rt(total_n * nbasis, df = 5) * sqrt(3 / 5), nrow = total_n, ncol = nbasis)
  }
  coef <- sweep(coef, 2, sd_pattern, FUN = "*")
  curves <- coef %*% t(basis_matrix)
  if (!is.null(mean_curve)) curves <- sweep(curves, 2, mean_curve, FUN = "+")
  curves
}

generate_ma1_curves <- function(total_n, basis_matrix, sd_pattern, kappa = 0.7, mean_curve = NULL) {
  nbasis <- ncol(basis_matrix)
  innov <- matrix(stats::rnorm((total_n + 1L) * nbasis), nrow = total_n + 1L, ncol = nbasis)
  innov <- sweep(innov, 2, sd_pattern, FUN = "*")
  coef <- innov[2:(total_n + 1L), , drop = FALSE] + kappa * innov[1:total_n, , drop = FALSE]
  curves <- coef %*% t(basis_matrix)
  if (!is.null(mean_curve)) curves <- sweep(curves, 2, mean_curve, FUN = "+")
  curves
}
