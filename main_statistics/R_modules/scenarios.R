# ==============================================================
# scenarios.R -- DGP generation and alternative changes
# ============================================================== 

generate_base_fd <- function(dgp_type, total_n) {
  baseInfo <- list(
    n = total_n,
    basisType = "bspline",
    muInfo = list(type = "sin", a = 1),
    factor = 1
  )
  if (dgp_type == "BB") {
    dataInfo <- c(baseInfo, list(nArgvals = 300))
    fd_data <- BB(dataInfo)
  } else if (dgp_type == "fIID") {
    dataInfo <- c(baseInfo, list(varType = "A", gaussian = TRUE))
    fd_data <- fIID(dataInfo)
  } else if (dgp_type == "fMA1") {
    dataInfo <- c(baseInfo, list(varType = "A", kappa = 0.7))
    fd_data <- fMA1(dataInfo)
  } else {
    stop(sprintf("Unsupported DGP type: %s", dgp_type))
  }
  if (!inherits(fd_data, "fd")) stop(sprintf("The object returned by %s is not of class 'fd'", dgp_type))
  fd_data
}

apply_change_to_fd <- function(fd_data,
                               scenario = c("null", "level_shift", "smooth_change", "abrupt_local_change", "sinusoidal_change"),
                               m,
                               T,
                               delta = 0,
                               s_star = 50L,
                               basis_k = localized_change_basis_index_default()) {
  scenario <- match.arg(scenario)
  total_n <- ncol(fd_data$coefs)
  if (scenario == "null" || isTRUE(all.equal(delta, 0))) return(fd_data)

  t_star <- m + as.integer(s_star)
  if (t_star > total_n) return(fd_data)
  change_idx <- t_star:total_n

  if (scenario == "level_shift") {
    fd_data$coefs[, change_idx] <- fd_data$coefs[, change_idx, drop = FALSE] + delta
    return(fd_data)
  }

  if (scenario == "smooth_change") {
    horizon_denom <- max(as.numeric(m) * as.numeric(T), 1)
    drift_weights <- delta * ((change_idx - t_star) / horizon_denom)
    fd_data$coefs[, change_idx] <- sweep(fd_data$coefs[, change_idx, drop = FALSE], 2, drift_weights, FUN = "+")
    return(fd_data)
  }

  if (scenario == "sinusoidal_change") {
    horizon_denom <- max(as.numeric(m) * as.numeric(T), 1)
    drift_weights <- delta * sin(pi * ((change_idx - t_star) / horizon_denom))
    fd_data$coefs[, change_idx] <- sweep(fd_data$coefs[, change_idx, drop = FALSE], 2, drift_weights, FUN = "+")
    return(fd_data)
  }

  if (scenario == "abrupt_local_change") {
    # The manuscript uses D = 21 spline basis functions overall and k = 5
    # as the default localized-break direction inside that basis.
    basis_dim <- streaming_curve_basis_dimension()
    basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = basis_dim)
    bk_coefs <- matrix(0, nrow = basis_dim, ncol = 1)
    basis_k <- min(max(1L, as.integer(basis_k)), basis_dim)
    bk_coefs[basis_k, 1] <- delta
    mu_fd <- fda::fd(coef = bk_coefs, basisobj = basis)
    for (i in change_idx) {
      fd_data$coefs[, i] <- fd_data$coefs[, i] + mu_fd$coefs[, 1]
    }
    return(fd_data)
  }

  fd_data
}
