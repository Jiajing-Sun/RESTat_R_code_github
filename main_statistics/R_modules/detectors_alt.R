# ==============================================================
# detectors_alt.R -- benchmark detector shapes after FPCA compression
# ============================================================== 

evaluate_alt_method <- function(method_row,
                                context_map,
                                q_used,
                                T_value,
                                alpha_levels,
                                cv_alt,
                                page_length_grid_size = 40L,
                                weighted_length_grid_size = 40L,
                                finite_eval_grid_size = 250L,
                                exact_page_scan = FALSE,
                                exact_weighted_scan = FALSE,
                                multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20))) {
  ctx <- context_map[[toupper(method_row$standardizer)]]
  proc <- ctx$process
  eval_idx <- build_eval_index_grid(proc$k_max, finite_eval_grid_size)
  U0 <- proc$U0_full[, 1:q_used, drop = FALSE]
  s <- proc$s
  s_full <- c(0, s)
  detector <- method_row$detector

  if (detector == "FullCUSUM") {
    D <- U0[eval_idx + 1L, , drop = FALSE]
    path <- quad_metric_rows(D, ctx) / g_gamma_sq(s[eval_idx], method_row$gamma)
    stat <- max(path)
    k_index <- eval_idx
  } else if (detector == "PageCUSUM") {
    lag_grid <- full_or_geometric_lag_grid(proc$k_max, page_length_grid_size, exact = exact_page_scan)
    best <- rep(-Inf, length(eval_idx))
    for (L in lag_grid) {
      hi <- eval_idx[eval_idx >= L]
      if (length(hi) == 0L) next
      pos <- match(hi, eval_idx)
      lo <- hi - L
      D <- U0[hi + 1L, , drop = FALSE] - U0[lo + 1L, , drop = FALSE]
      lens <- s_full[hi + 1L] - s_full[lo + 1L]
      cur <- quad_metric_rows(D, ctx) / g_gamma_sq(lens, method_row$gamma)
      best[pos] <- pmax(best[pos], cur)
    }
    path <- best
    stat <- max(path)
    k_index <- eval_idx
  } else if (detector == "WeightedCUSUM") {
    lag_grid <- full_or_geometric_lag_grid(proc$k_max, weighted_length_grid_size, exact = exact_weighted_scan)
    best <- rep(-Inf, length(eval_idx))
    for (L in lag_grid) {
      hi <- eval_idx[eval_idx >= L]
      if (length(hi) == 0L) next
      pos <- match(hi, eval_idx)
      lo <- hi - L
      D <- U0[hi + 1L, , drop = FALSE] - U0[lo + 1L, , drop = FALSE]
      lens <- s_full[hi + 1L] - s_full[lo + 1L]
      base <- quad_metric_rows(D, ctx) / g_gamma_sq(lens, method_row$gamma)
      cur <- base * make_weighted_cusum_omega(lens, method_row$omega_name)
      best[pos] <- pmax(best[pos], cur)
    }
    path <- best
    stat <- max(path)
    k_index <- eval_idx
  } else if (detector == "MOSUM") {
    h <- as.numeric(method_row$bandwidth_h)
    G <- as.integer(floor(proc$m * h))
    valid_hi <- eval_idx[eval_idx >= 2L * G]
    D <- U0[valid_hi + 1L, , drop = FALSE] - 2 * U0[valid_hi - G + 1L, , drop = FALSE] + U0[valid_hi - 2L * G + 1L, , drop = FALSE]
    path <- quad_metric_rows(D, ctx) / h
    stat <- max(path)
    k_index <- valid_hi
  } else if (detector == "MultiscaleMOSUM") {
    h_vec <- multiscale_h_sets[[method_row$hset_name]]
    if (is.null(h_vec)) stop("Unknown multiscale h-set: ", method_row$hset_name)
    best <- rep(-Inf, length(eval_idx))
    for (h in h_vec) {
      G <- as.integer(floor(proc$m * h))
      valid_hi <- eval_idx[eval_idx >= 2L * G]
      if (length(valid_hi) == 0L) next
      pos <- match(valid_hi, eval_idx)
      D <- U0[valid_hi + 1L, , drop = FALSE] - 2 * U0[valid_hi - G + 1L, , drop = FALSE] + U0[valid_hi - 2L * G + 1L, , drop = FALSE]
      cur <- quad_metric_rows(D, ctx) / h
      cur <- cur * make_multiscale_weight(h, method_row$scale_weight_name)
      best[pos] <- pmax(best[pos], cur)
    }
    path <- best
    stat <- max(path)
    k_index <- eval_idx
  } else {
    stop("Unknown benchmark detector: ", detector)
  }

  out <- vector("list", length(alpha_levels))
  for (i in seq_along(alpha_levels)) {
    a <- alpha_levels[i]
    cv <- lookup_alt_critical_value(cv_alt, standardizer = method_row$standardizer, detector = method_row$detector,
                                    T = T_value, q = q_used, alpha = a, gamma = method_row$gamma,
                                    bandwidth_h = method_row$bandwidth_h, omega_name = method_row$omega_name,
                                    hset_name = method_row$hset_name, scale_weight_name = method_row$scale_weight_name)
    out[[i]] <- data.frame(
      family = method_row$family,
      standardizer = method_row$standardizer,
      detector = method_row$detector,
      type = method_row$type,
      gamma = method_row$gamma,
      weight_name = NA_character_,
      bandwidth_h = method_row$bandwidth_h,
      omega_name = ifelse(method_row$omega_name == "", NA, method_row$omega_name),
      hset_name = ifelse(method_row$hset_name == "", NA, method_row$hset_name),
      scale_weight_name = ifelse(method_row$scale_weight_name == "", NA, method_row$scale_weight_name),
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
