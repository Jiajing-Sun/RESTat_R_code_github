# ==============================================================
# contamination.R -- Phase-I contamination mechanisms for the
# streaming-curve simulations
# ============================================================== 

contamination_shape_values <- function(t_grid,
                                       basis_matrix = NULL,
                                       scenario_ref = c("level_shift", "smooth_change", "abrupt_local_change"),
                                       basis_k = 5L) {
  scenario_ref <- match.arg(scenario_ref)
  if (scenario_ref %in% c("level_shift", "smooth_change")) {
    return(rep(1, length(t_grid)))
  }
  local_bump_shape(t_grid, basis_matrix = basis_matrix, basis_k = basis_k)
}

apply_training_contamination_to_curves <- function(curves,
                                                   t_grid,
                                                   basis_matrix = NULL,
                                                   m,
                                                   contam_type = c("none", "late_break", "drift"),
                                                   train_b = 0,
                                                   train_break_frac = 0.8,
                                                   train_drift_start_frac = 0.0,
                                                   contam_shape = c("phase2_matched", "level_shift", "abrupt_local_change"),
                                                   phase2_scenario = c("level_shift", "smooth_change", "abrupt_local_change"),
                                                   basis_k = 5L) {
  contam_type <- match.arg(contam_type)
  curves <- safe_matrix(curves)
  if (contam_type == "none" || isTRUE(all.equal(as.numeric(train_b), 0))) {
    return(list(curves = curves, t_train_star = NA_integer_, contam_shape_name = "none"))
  }

  contam_shape <- match.arg(contam_shape)
  phase2_scenario <- match.arg(phase2_scenario)
  shape_scenario <- if (contam_shape == "phase2_matched") phase2_scenario else contam_shape
  shift_shape <- contamination_shape_values(t_grid, basis_matrix = basis_matrix,
                                            scenario_ref = shape_scenario,
                                            basis_k = basis_k)

  total_n <- nrow(curves)
  if (m < 2L) stop("m must be at least 2 for training contamination designs")
  if (m > total_n) stop("m exceeds the number of curves")

  if (contam_type == "late_break") {
    t_star <- as.integer(floor(train_break_frac * m))
    t_star <- max(1L, min(t_star, m - 1L))
    idx <- (t_star + 1L):total_n
    curves[idx, ] <- sweep(curves[idx, , drop = FALSE], 2, as.numeric(train_b) * shift_shape, FUN = "+")
    return(list(curves = curves, t_train_star = t_star, contam_shape_name = shape_scenario))
  }

  if (contam_type == "drift") {
    t0 <- as.integer(floor(train_drift_start_frac * m))
    t0 <- max(0L, min(t0, m - 1L))
    t_start <- t0 + 1L
    idx_train <- t_start:m
    if (length(idx_train) == 1L) {
      drift_weights <- as.numeric(train_b)
    } else {
      drift_weights <- as.numeric(train_b) * (idx_train - t_start) / (m - t_start)
    }
    for (i in seq_along(idx_train)) {
      curves[idx_train[i], ] <- curves[idx_train[i], ] + drift_weights[i] * shift_shape
    }
    if (m < total_n) {
      idx_post <- (m + 1L):total_n
      curves[idx_post, ] <- sweep(curves[idx_post, , drop = FALSE], 2, as.numeric(train_b) * shift_shape, FUN = "+")
    }
    return(list(curves = curves, t_train_star = t_start, contam_shape_name = shape_scenario))
  }

  stop("Unknown contamination type: ", contam_type)
}
