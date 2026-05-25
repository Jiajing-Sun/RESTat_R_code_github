# ==============================================================
# genData.R -- minimal DGP generators used by the simulation code
# All fda functions are namespace-qualified so PSOCK workers do not
# depend on library(fda) being attached.
# ==============================================================

make_basis_obj <- function(basis_type = c("bspline", "fourier"), nbasis = NULL, n_breaks = 301L) {
  basis_type <- match.arg(basis_type)
  if (is.null(nbasis)) nbasis <- streaming_curve_basis_dimension()
  if (basis_type == "bspline") {
    fda::create.bspline.basis(rangeval = c(0, 1), nbasis = as.integer(nbasis))
  } else {
    fda::create.fourier.basis(rangeval = c(0, 1), nbasis = as.integer(nbasis))
  }
}

make_mu_values <- function(t_grid, muInfo) {
  mu_type <- muInfo$type %||% "sin"
  mu_a <- as.numeric(muInfo$a %||% 0)
  if (identical(mu_a, 0)) return(rep(0, length(t_grid)))
  if (mu_type == "sin") return(mu_a * sin(2 * pi * t_grid))
  if (mu_type == "horv") return(mu_a * t_grid * (1 - t_grid))
  stop("Unsupported muInfo$type: ", mu_type)
}

make_mu_fd <- function(basis, muInfo, n_grid = 301L) {
  t_grid <- seq(0, 1, length.out = as.integer(n_grid))
  mu_vals <- make_mu_values(t_grid, muInfo)
  fda::Data2fd(argvals = t_grid, y = mu_vals, basisobj = basis, lambda = 0)
}

BB <- function(dataInfo) {
  n <- as.integer(dataInfo$n)
  N <- as.integer(dataInfo$nArgvals %||% 300L)
  muInfo <- dataInfo$muInfo %||% list(type = "sin", a = 0)
  basis_type <- dataInfo$basisType %||% "bspline"

  eps <- matrix(stats::rnorm(N * n), nrow = N, ncol = n)
  values <- matrix(0, nrow = N + 1L, ncol = n)
  t_grid <- seq(0, 1, length.out = N + 1L)
  t_no0 <- t_grid[-1L]
  scale <- 1 / sqrt(N)
  for (i in seq_len(n)) {
    cumulative <- scale * cumsum(eps[, i])
    values[-1L, i] <- cumulative - t_no0 * cumulative[N]
  }
  mu_vals <- make_mu_values(t_grid, muInfo)
  values <- sweep(values, 1L, mu_vals, FUN = "+")

  basis <- if (basis_type == "bspline") {
    fda::create.bspline.basis(rangeval = c(0, 1), norder = 2L, breaks = t_grid)
  } else {
    fda::create.fourier.basis(rangeval = c(0, 1), nbasis = N + 1L)
  }

  fda::Data2fd(argvals = t_grid, y = values, basisobj = basis, lambda = 0)
}

fIID <- function(dataInfo) {
  n <- as.integer(dataInfo$n)
  muInfo <- dataInfo$muInfo %||% list(type = "sin", a = 0)
  factor <- as.numeric(dataInfo$factor %||% 1)
  gaussian <- if (is.null(dataInfo$gaussian)) TRUE else isTRUE(dataInfo$gaussian)
  var_type <- dataInfo$varType %||% "A"
  basis_type <- dataInfo$basisType %||% "bspline"

  basis <- make_basis_obj(basis_type, nbasis = streaming_curve_basis_dimension())
  D <- basis$nbasis
  Sigma <- switch(var_type,
                  A = 1 / seq_len(D),
                  B = 1.2^(-seq_len(D)),
                  stop("Unsupported varType: ", var_type))

  coef <- matrix(0, nrow = D, ncol = n)
  for (i in seq_len(n)) {
    if (gaussian) {
      coef[, i] <- stats::rnorm(D, mean = 0, sd = Sigma)
    } else {
      coef[, i] <- stats::rt(D, df = 5) * sqrt(3 / 5) * Sigma
    }
  }

  fdata <- fda::fd(coef = factor * coef, basisobj = basis)
  if (!is.null(muInfo) && as.numeric(muInfo$a %||% 0) != 0) {
    mu_fd <- make_mu_fd(basis, muInfo)
    fdata$coefs <- sweep(fdata$coefs, 1L, mu_fd$coefs[, 1L], FUN = "+")
  }
  fdata
}

createMA <- function(n, basis, Theta, Sigma) {
  n <- as.integer(n)
  D <- basis$nbasis
  zlag0 <- matrix(0, nrow = D, ncol = n + 1L)
  for (i in seq_len(n + 1L)) zlag0[, i] <- stats::rnorm(D, mean = 0, sd = Sigma)
  zlag1 <- matrix(0, nrow = D, ncol = n)
  if (n >= 2L) {
    for (i in 2:n) zlag1[, i] <- Theta %*% zlag0[, i - 1L]
  }
  coef <- zlag1 + zlag0[, 2:(n + 1L), drop = FALSE]
  fda::fd(coef = coef, basisobj = basis)
}

fMA1 <- function(dataInfo) {
  n <- as.integer(dataInfo$n)
  muInfo <- dataInfo$muInfo %||% list(type = "sin", a = 0)
  factor <- as.numeric(dataInfo$factor %||% 1)
  var_type <- dataInfo$varType %||% "A"
  basis_type <- dataInfo$basisType %||% "bspline"
  kappa <- as.numeric(dataInfo$kappa %||% 0.7)

  basis <- make_basis_obj(basis_type, nbasis = streaming_curve_basis_dimension())
  D <- basis$nbasis
  Sigma <- switch(var_type,
                  A = 1 / seq_len(D),
                  B = 1.2^(-seq_len(D)),
                  stop("Unsupported varType: ", var_type))

  Psi <- dataInfo$Psi
  if (is.null(Psi)) {
    Psi <- matrix(0, nrow = D, ncol = D)
    for (i in seq_len(D)) {
      for (j in seq_len(D)) {
        Psi[i, j] <- stats::rnorm(1L, mean = 0, sd = Sigma[i] * Sigma[j])
      }
    }
    ev1 <- max(Re(eigen(Psi %*% t(Psi), symmetric = TRUE, only.values = TRUE)$values))
    Psi <- Psi / sqrt(max(ev1, 1e-8))
  }

  Theta <- kappa * Psi
  fdata <- createMA(n = n, basis = basis, Theta = Theta, Sigma = Sigma)
  fdata$coefs <- fdata$coefs * factor

  if (!is.null(muInfo) && as.numeric(muInfo$a %||% 0) != 0) {
    mu_fd <- make_mu_fd(basis, muInfo)
    fdata$coefs <- sweep(fdata$coefs, 1L, mu_fd$coefs[, 1L], FUN = "+")
  }
  fdata
}
