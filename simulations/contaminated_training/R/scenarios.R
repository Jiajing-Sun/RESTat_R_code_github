# ==============================================================
# scenarios.R -- DGP generation and structural changes on a dense
# grid representation of the curves
# ============================================================== 

generate_base_curves <- function(dgp_type,
                                 total_n,
                                 n_grid = 301L,
                                 basis_type = "bspline",
                                 nbasis = 21L,
                                 muInfo = list(type = "sin", a = 1),
                                 varType = "A",
                                 gaussian = TRUE,
                                 kappa = 0.7) {
  t_grid <- seq(0, 1, length.out = as.integer(n_grid))
  basis_matrix <- make_basis_matrix(t_grid, nbasis = nbasis, basis_type = basis_type)
  mean_curve <- make_mu_values(t_grid, muInfo)
  sd_pattern <- switch(varType,
                       A = 1 / seq_len(nbasis),
                       B = 1.2^(-seq_len(nbasis)),
                       stop("Unsupported varType: ", varType))

  values <- if (dgp_type == "BB") {
    generate_bb_curves(total_n = total_n, t_grid = t_grid, mean_curve = mean_curve)
  } else if (dgp_type == "fIID") {
    generate_iid_curves(total_n = total_n, basis_matrix = basis_matrix, sd_pattern = sd_pattern,
                        mean_curve = mean_curve, gaussian = gaussian)
  } else if (dgp_type == "fMA1") {
    generate_ma1_curves(total_n = total_n, basis_matrix = basis_matrix, sd_pattern = sd_pattern,
                        kappa = kappa, mean_curve = mean_curve)
  } else {
    stop(sprintf("Unsupported DGP type: %s", dgp_type))
  }

  list(values = safe_matrix(values), t_grid = t_grid, basis_matrix = basis_matrix)
}

local_bump_shape <- function(t_grid, basis_matrix = NULL, basis_k = 5L) {
  if (!is.null(basis_matrix) && ncol(basis_matrix) >= 1L) {
    k <- min(max(1L, as.integer(basis_k)), ncol(basis_matrix))
    shp <- basis_matrix[, k]
    mx <- max(abs(shp))
    if (mx > 0) return(shp / mx)
  }
  center <- 0.5
  width <- 0.10
  shp <- exp(-0.5 * ((t_grid - center) / width)^2)
  shp / max(shp)
}

apply_change_to_curves <- function(curves,
                                   t_grid,
                                   basis_matrix = NULL,
                                   scenario = c("null", "level_shift", "smooth_change", "abrupt_local_change"),
                                   m,
                                   T,
                                   delta = 0,
                                   s_star = 50,
                                   basis_k = 5L) {
  scenario <- match.arg(scenario)
  curves <- safe_matrix(curves)
  total_n <- nrow(curves)
  if (scenario == "null" || isTRUE(all.equal(delta, 0))) return(curves)

  t_star <- m + as.integer(s_star)
  if (t_star > total_n) return(curves)
  change_idx <- t_star:total_n

  if (scenario == "level_shift") {
    shift_shape <- rep(1, ncol(curves))
    curves[change_idx, ] <- sweep(curves[change_idx, , drop = FALSE], 2, delta * shift_shape, FUN = "+")
    return(curves)
  }

  if (scenario == "smooth_change") {
    shift_shape <- rep(1, ncol(curves))
    num_change <- length(change_idx)
    drift_weights <- if (delta == 0) rep(0, num_change) else seq(0, delta, length.out = num_change)
    for (i in seq_along(change_idx)) {
      curves[change_idx[i], ] <- curves[change_idx[i], ] + drift_weights[i] * shift_shape
    }
    return(curves)
  }

  if (scenario == "abrupt_local_change") {
    bump <- local_bump_shape(t_grid, basis_matrix = basis_matrix, basis_k = basis_k)
    curves[change_idx, ] <- sweep(curves[change_idx, , drop = FALSE], 2, delta * bump, FUN = "+")
    return(curves)
  }

  curves
}
