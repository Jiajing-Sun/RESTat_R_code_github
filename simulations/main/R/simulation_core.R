
# ==============================================================
# simulation_core.R -- one-replicate simulation engine (robust)
# ============================================================== 

fill_or_create_column <- function(df, column, value) {
  n <- nrow(df)
  if (!column %in% names(df)) {
    df[[column]] <- rep(value, length.out = n)
    return(df)
  }
  idx <- is.na(df[[column]])
  if (is.character(df[[column]])) idx <- idx | !nzchar(df[[column]])
  if (any(idx)) df[[column]][idx] <- rep(value, length.out = sum(idx))
  df
}

make_method_error_rows <- function(method_row, alpha_levels, rep_id, sim_par, err_msg) {
  out <- vector("list", length(alpha_levels))
  for (i in seq_along(alpha_levels)) {
    a <- as.numeric(alpha_levels[i])
    out[[i]] <- data.frame(
      rep = rep_id,
      dgp_type = sim_par$dgp_type,
      scenario = sim_par$scenario,
      m = sim_par$m,
      T = canonical_T_scalar(sim_par$T_value),
      delta = sim_par$delta,
      s_star = sim_par$s_star,
      q_used = NA_integer_,
      fpca_mode = if (is.na(sim_par$fixed_q)) "FVE" else "fixed_q",
      fve_threshold = sim_par$fve_threshold,
      family = method_row$family %||% NA_character_,
      standardizer = method_row$standardizer %||% NA_character_,
      detector = method_row$detector %||% NA_character_,
      type = method_row$type %||% NA_character_,
      gamma = if ("gamma" %in% names(method_row)) as.numeric(method_row$gamma) else NA_real_,
      weight_name = if ("weight_name" %in% names(method_row) && nzchar(method_row$weight_name)) method_row$weight_name else NA_character_,
      bandwidth_h = if ("bandwidth_h" %in% names(method_row)) as.numeric(method_row$bandwidth_h) else NA_real_,
      omega_name = if ("omega_name" %in% names(method_row) && nzchar(method_row$omega_name)) method_row$omega_name else NA_character_,
      hset_name = if ("hset_name" %in% names(method_row) && nzchar(method_row$hset_name)) method_row$hset_name else NA_character_,
      scale_weight_name = if ("scale_weight_name" %in% names(method_row) && nzchar(method_row$scale_weight_name)) method_row$scale_weight_name else NA_character_,
      method_group = method_row$method_group %||% NA_character_,
      method_label = method_row$method_label %||% NA_character_,
      method_id = method_row$method_id %||% NA_character_,
      alpha = a,
      statistic = NA_real_,
      critical_value = NA_real_,
      reject = NA,
      first_rejection = NA_integer_,
      error_flag = TRUE,
      error_message = as.character(err_msg),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

make_global_error_rows <- function(rep_id, sim_par, alpha_levels, err_msg) {
  out <- vector("list", length(alpha_levels))
  for (i in seq_along(alpha_levels)) {
    a <- as.numeric(alpha_levels[i])
    out[[i]] <- data.frame(
      rep = rep_id,
      dgp_type = sim_par$dgp_type,
      scenario = sim_par$scenario,
      m = sim_par$m,
      T = canonical_T_scalar(sim_par$T_value),
      delta = sim_par$delta,
      s_star = sim_par$s_star,
      q_used = NA_integer_,
      fpca_mode = if (is.na(sim_par$fixed_q)) "FVE" else "fixed_q",
      fve_threshold = sim_par$fve_threshold,
      family = NA_character_,
      standardizer = NA_character_,
      detector = NA_character_,
      type = NA_character_,
      gamma = NA_real_,
      weight_name = NA_character_,
      bandwidth_h = NA_real_,
      omega_name = NA_character_,
      hset_name = NA_character_,
      scale_weight_name = NA_character_,
      method_group = "GlobalError",
      method_label = "GlobalError",
      method_id = NA_character_,
      alpha = a,
      statistic = NA_real_,
      critical_value = NA_real_,
      reject = NA,
      first_rejection = NA_integer_,
      error_flag = TRUE,
      error_message = as.character(err_msg),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

safe_eval_main_method <- function(method_row, context_map, q_used, T_value, alpha_levels, cv_main, rep_id, sim_par) {
  tryCatch(
    evaluate_main_method(method_row, context_map, q_used, T_value = T_value, alpha_levels = alpha_levels, cv_main = cv_main),
    error = function(e) make_method_error_rows(method_row, alpha_levels, rep_id, sim_par, conditionMessage(e))
  )
}

safe_eval_alt_method <- function(method_row, context_map, q_used, T_value, alpha_levels, cv_alt,
                                 page_length_grid_size, weighted_length_grid_size, finite_eval_grid_size,
                                 exact_page_scan, exact_weighted_scan, multiscale_h_sets,
                                 rep_id, sim_par) {
  tryCatch(
    evaluate_alt_method(method_row, context_map, q_used, T_value = T_value, alpha_levels = alpha_levels,
                        cv_alt = cv_alt,
                        page_length_grid_size = page_length_grid_size,
                        weighted_length_grid_size = weighted_length_grid_size,
                        finite_eval_grid_size = finite_eval_grid_size,
                        exact_page_scan = exact_page_scan,
                        exact_weighted_scan = exact_weighted_scan,
                        multiscale_h_sets = multiscale_h_sets),
    error = function(e) make_method_error_rows(method_row, alpha_levels, rep_id, sim_par, conditionMessage(e))
  )
}

simulate_one_replicate <- function(rep_id,
                                   sim_par,
                                   cv_main,
                                   cv_alt,
                                   main_catalog,
                                   alt_catalog) {
  tryCatch({
    total_n <- if (is.finite(sim_par$T_value)) {
      sim_par$m + as.integer(round(sim_par$m * sim_par$T_value))
    } else {
      stop("The current simulation drivers focus on finite-horizon designs. Use finite T in the parameter grid.")
    }

    fd_data <- generate_base_fd(sim_par$dgp_type, total_n = total_n)
    fd_data <- apply_change_to_fd(fd_data,
                                  scenario = sim_par$scenario,
                                  m = sim_par$m,
                                  T = sim_par$T_value,
                                  delta = sim_par$delta,
                                  s_star = sim_par$s_star,
                                  basis_k = sim_par$basis_k)

    fpca <- fpca_project_scores(fd_data,
                                m = sim_par$m,
                                q_cap = sim_par$q_cap,
                                fve_threshold = sim_par$fve_threshold,
                                fixed_q = sim_par$fixed_q,
                                n_grid = sim_par$n_grid)

    q_used <- fpca$q_used
    whitened <- whiten_scores(fpca$scores_train, fpca$scores_monitor, ridge = sim_par$ridge)

    ctx_raw_hac <- prepare_standardizer_context(fpca$scores_train, fpca$scores_monitor, "HAC",
                                                ridge = sim_par$ridge,
                                                range_floor = sim_par$range_floor,
                                                hac_bandwidth = sim_par$hac_bandwidth)
    ctx_raw_ss <- prepare_standardizer_context(fpca$scores_train, fpca$scores_monitor, "SSMS",
                                               ridge = sim_par$ridge,
                                               range_floor = sim_par$range_floor,
                                               hac_bandwidth = sim_par$hac_bandwidth)
    ctx_pw_rs <- prepare_standardizer_context(whitened$scores_train, whitened$scores_monitor, "RSMS",
                                              ridge = sim_par$ridge,
                                              range_floor = sim_par$range_floor,
                                              hac_bandwidth = sim_par$hac_bandwidth)

    context_map <- list(HAC = ctx_raw_hac, SSMS = ctx_raw_ss, RSMS = ctx_pw_rs)

    out_list <- list()
    idx <- 1L

    if (nrow(main_catalog) > 0L) {
      for (j in seq_len(nrow(main_catalog))) {
        out_list[[idx]] <- safe_eval_main_method(main_catalog[j, , drop = FALSE], context_map, q_used,
                                                 T_value = sim_par$T_value, alpha_levels = sim_par$alpha_levels,
                                                 cv_main = cv_main, rep_id = rep_id, sim_par = sim_par)
        idx <- idx + 1L
      }
    }

    if (nrow(alt_catalog) > 0L) {
      for (j in seq_len(nrow(alt_catalog))) {
        out_list[[idx]] <- safe_eval_alt_method(alt_catalog[j, , drop = FALSE], context_map, q_used,
                                                T_value = sim_par$T_value, alpha_levels = sim_par$alpha_levels,
                                                cv_alt = cv_alt,
                                                page_length_grid_size = sim_par$page_length_grid_size,
                                                weighted_length_grid_size = sim_par$weighted_length_grid_size,
                                                finite_eval_grid_size = sim_par$finite_eval_grid_size,
                                                exact_page_scan = sim_par$exact_page_scan,
                                                exact_weighted_scan = sim_par$exact_weighted_scan,
                                                multiscale_h_sets = sim_par$multiscale_h_sets,
                                                rep_id = rep_id, sim_par = sim_par)
        idx <- idx + 1L
      }
    }

    out <- do.call(rbind, out_list)
    out <- fill_or_create_column(out, "rep", rep_id)
    out <- fill_or_create_column(out, "dgp_type", sim_par$dgp_type)
    out <- fill_or_create_column(out, "scenario", sim_par$scenario)
    out <- fill_or_create_column(out, "m", sim_par$m)
    out <- fill_or_create_column(out, "T", canonical_T_scalar(sim_par$T_value))
    out <- fill_or_create_column(out, "delta", sim_par$delta)
    out <- fill_or_create_column(out, "s_star", sim_par$s_star)
    out <- fill_or_create_column(out, "q_used", q_used)
    out <- fill_or_create_column(out, "fpca_mode", if (is.na(sim_par$fixed_q)) "FVE" else "fixed_q")
    out <- fill_or_create_column(out, "fve_threshold", sim_par$fve_threshold)
    out <- fill_or_create_column(out, "error_flag", FALSE)
    out <- fill_or_create_column(out, "error_message", "")
    out
  }, error = function(e) {
    make_global_error_rows(rep_id, sim_par, sim_par$alpha_levels, conditionMessage(e))
  })
}
