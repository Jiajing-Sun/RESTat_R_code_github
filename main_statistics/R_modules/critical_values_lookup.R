# ==============================================================
# critical_values_lookup.R -- load and query the latest critical values
# ============================================================== 

file_has_nul_bytes <- function(path, chunk_size = 65536L) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  repeat {
    buf <- readBin(con, what = "raw", n = chunk_size)
    if (length(buf) == 0L) return(FALSE)
    if (any(buf == as.raw(0))) return(TRUE)
  }
}

read_csv_strict <- function(path) {
  if (file_has_nul_bytes(path)) {
    stop("Detected NUL bytes in CSV file: ", path)
  }

  warnings_seen <- character()
  out <- withCallingHandlers(
    read.csv(path, stringsAsFactors = FALSE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (length(warnings_seen) > 0L) {
    stop(
      "Warnings were raised while reading CSV file ", path, ": ",
      paste(unique(warnings_seen), collapse = " | ")
    )
  }

  out
}

candidate_main_cv_paths <- function(root) {
  list(
    base = unique(c(
      file.path(root, "critical_values", "critical_values_all.csv"),
      file.path(root, "streaming_curve_cv_codes_windowsfix_v3_hac", "outputs", "critical_values_all.csv")
    )),
    weights = unique(c(
      file.path(root, "critical_values", "critical_values_all_weights.csv"),
      file.path(root, "streaming_curve_cv_codes_windowsfix_v3_hac", "outputs", "critical_values_all_weights.csv")
    ))
  )
}

candidate_alt_cv_paths <- function(root) {
  unique(c(
    file.path(root, "critical_values", "critical_values_alt_detectors.csv"),
    file.path(root, "critical_values_alt_detectors.csv"),
    file.path(root, "alt_detector_cv_codes_fastfix_v2", "outputs", "critical_values_alt_detectors.csv")
  ))
}

read_first_valid_csv <- function(paths, label) {
  errs <- character()
  for (path in unique(paths)) {
    if (!file.exists(path)) next
    out <- tryCatch(read_csv_strict(path), error = function(e) e)
    if (!inherits(out, "error")) {
      return(list(data = out, path = path, errors = errs))
    }
    errs <- c(errs, sprintf("%s :: %s", normalize_path2(path, mustWork = FALSE), conditionMessage(out)))
  }

  stop(
    "Could not load ", label, ". Tried these candidates:\n",
    paste(unique(c(paths, errs)), collapse = "\n")
  )
}

normalize_main_critical_values <- function(base, weights) {
  base$T_label <- canonical_T_vec(base$T)
  weights$T_label <- canonical_T_vec(weights$T)
  base$stat_upper <- toupper(base$stat)
  base$type_upper <- toupper(base$type)
  weights$stat_upper <- toupper(weights$stat)
  weights$type_upper <- toupper(weights$type)
  weights$weight_norm <- vapply(weights$weight_name, normalize_weight_name, character(1L))
  list(base = base, weights = weights)
}

normalize_alt_critical_values <- function(df) {
  df$T_label <- canonical_T_vec(df$T)
  df$standardizer_upper <- toupper(df$standardizer)
  df$detector_upper <- toupper(df$detector)
  if ("omega_name" %in% names(df)) {
    df$omega_norm <- ifelse(is.na(df$omega_name) | df$omega_name == "", "", vapply(ifelse(is.na(df$omega_name), "", df$omega_name), function(z) if (z == "") "" else normalize_omega_name(z), character(1L)))
  } else {
    df$omega_norm <- ""
  }
  if ("scale_weight_name" %in% names(df)) {
    df$scale_weight_norm <- ifelse(is.na(df$scale_weight_name) | df$scale_weight_name == "", "", vapply(ifelse(is.na(df$scale_weight_name), "", df$scale_weight_name), function(z) if (z == "") "" else normalize_scale_weight_name(z), character(1L)))
  } else {
    df$scale_weight_norm <- ""
  }
  df
}

load_main_critical_values <- function(root) {
  paths <- candidate_main_cv_paths(root)
  base_result <- read_first_valid_csv(paths$base, label = "main critical-value file")
  weight_result <- read_first_valid_csv(paths$weights, label = "weighted critical-value file")
  if (!identical(base_result$path, paths$base[1L])) {
    message("Using fallback main critical-value file: ", normalize_path2(base_result$path, mustWork = FALSE))
  }
  if (!identical(weight_result$path, paths$weights[1L])) {
    message("Using fallback weighted critical-value file: ", normalize_path2(weight_result$path, mustWork = FALSE))
  }
  base <- base_result$data
  weights <- weight_result$data
  normalize_main_critical_values(base, weights)
}

load_alt_critical_values <- function(root) {
  result <- read_first_valid_csv(candidate_alt_cv_paths(root), label = "alternative-detector critical-value file")
  if (!identical(result$path, candidate_alt_cv_paths(root)[1L])) {
    message("Using fallback benchmark critical-value file: ", normalize_path2(result$path, mustWork = FALSE))
  }
  df <- result$data
  normalize_alt_critical_values(df)
}

lookup_main_critical_value <- function(cv_main, standardizer, type, T, q, alpha, gamma = NA_real_, weight_name = "U") {
  T_label <- canonical_T_scalar(T)
  q <- as.integer(q)
  alpha <- as.numeric(alpha)
  stat_upper <- toupper(standardizer)
  type_upper <- toupper(type)

  if (type_upper == "KS") {
    df <- cv_main$base
    keep <- df$stat_upper == stat_upper &
      df$type_upper == "KS" &
      df$T_label == T_label &
      df$q == q &
      abs(df$alpha - alpha) < 1e-12 &
      abs(df$gamma - as.numeric(gamma)) < 1e-12
    sub <- df[keep, , drop = FALSE]
  } else {
    wnorm <- normalize_weight_name(weight_name)
    if (wnorm == "U") {
      df <- cv_main$base
      keep <- df$stat_upper == stat_upper &
        df$type_upper == "CVM" &
        df$T_label == T_label &
        df$q == q &
        abs(df$alpha - alpha) < 1e-12
      sub <- df[keep, , drop = FALSE]
    } else {
      df <- cv_main$weights
      keep <- df$stat_upper == stat_upper &
        df$type_upper == "CVM" &
        df$T_label == T_label &
        df$q == q &
        abs(df$alpha - alpha) < 1e-12 &
        df$weight_norm == wnorm
      sub <- df[keep, , drop = FALSE]
    }
  }

  if (nrow(sub) != 1L) {
    stop(sprintf("Main critical value not uniquely found for (%s,%s,T=%s,q=%s,alpha=%s,gamma=%s,weight=%s). Matches=%d",
                 standardizer, type, T_label, q, alpha, gamma, weight_name, nrow(sub)))
  }
  as.numeric(sub$critical_value[1L])
}

collect_main_critical_value_coverage_issues <- function(cv_main, main_catalog, T_grid, q_grid, alpha_levels) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(main_catalog))) {
    row <- main_catalog[i, , drop = FALSE]
    for (T_value in T_grid) {
      for (q in q_grid) {
        for (alpha in alpha_levels) {
          err <- tryCatch({
            lookup_main_critical_value(
              cv_main,
              standardizer = row$standardizer,
              type = row$type,
              T = T_value,
              q = q,
              alpha = alpha,
              gamma = row$gamma,
              weight_name = row$weight_name
            )
            NULL
          }, error = function(e) conditionMessage(e))
          if (!is.null(err)) {
            rows[[idx]] <- data.frame(
              family = row$family,
              standardizer = row$standardizer,
              detector = row$detector,
              type = row$type,
              gamma = row$gamma,
              weight_name = row$weight_name,
              T = canonical_T_scalar(T_value),
              q = as.integer(q),
              alpha = as.numeric(alpha),
              issue = err,
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

collect_alt_critical_value_coverage_issues <- function(cv_alt, alt_catalog, T_grid, q_grid, alpha_levels) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(alt_catalog))) {
    row <- alt_catalog[i, , drop = FALSE]
    for (T_value in T_grid) {
      for (q in q_grid) {
        for (alpha in alpha_levels) {
          err <- tryCatch({
            lookup_alt_critical_value(
              cv_alt,
              standardizer = row$standardizer,
              detector = row$detector,
              T = T_value,
              q = q,
              alpha = alpha,
              gamma = row$gamma,
              bandwidth_h = row$bandwidth_h,
              omega_name = row$omega_name,
              hset_name = row$hset_name,
              scale_weight_name = row$scale_weight_name
            )
            NULL
          }, error = function(e) conditionMessage(e))
          if (!is.null(err)) {
            rows[[idx]] <- data.frame(
              family = row$family,
              standardizer = row$standardizer,
              detector = row$detector,
              type = row$type,
              gamma = row$gamma,
              bandwidth_h = row$bandwidth_h,
              omega_name = row$omega_name,
              hset_name = row$hset_name,
              scale_weight_name = row$scale_weight_name,
              T = canonical_T_scalar(T_value),
              q = as.integer(q),
              alpha = as.numeric(alpha),
              issue = err,
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

validate_main_critical_value_coverage <- function(cv_main, main_catalog, T_grid, q_grid, alpha_levels) {
  issues <- collect_main_critical_value_coverage_issues(cv_main, main_catalog, T_grid, q_grid, alpha_levels)
  if (nrow(issues) > 0L) {
    preview <- utils::capture.output(print(utils::head(issues, 10L)))
    stop(
      "Main critical-value coverage failed for the requested grid.\n",
      paste(preview, collapse = "\n")
    )
  }
  invisible(TRUE)
}

validate_alt_critical_value_coverage <- function(cv_alt, alt_catalog, T_grid, q_grid, alpha_levels) {
  issues <- collect_alt_critical_value_coverage_issues(cv_alt, alt_catalog, T_grid, q_grid, alpha_levels)
  if (nrow(issues) > 0L) {
    preview <- utils::capture.output(print(utils::head(issues, 10L)))
    stop(
      "Benchmark-detector critical-value coverage failed for the requested grid.\n",
      paste(preview, collapse = "\n")
    )
  }
  invisible(TRUE)
}

lookup_alt_critical_value <- function(cv_alt, standardizer, detector, T, q, alpha,
                                      gamma = NA_real_, bandwidth_h = NA_real_, omega_name = "",
                                      hset_name = "", scale_weight_name = "") {
  T_label <- canonical_T_scalar(T)
  q <- as.integer(q)
  alpha <- as.numeric(alpha)
  std_u <- toupper(standardizer)
  det_u <- toupper(detector)
  omega_norm <- if (!is.null(omega_name) && nzchar(omega_name)) normalize_omega_name(omega_name) else ""
  scale_norm <- if (!is.null(scale_weight_name) && nzchar(scale_weight_name)) normalize_scale_weight_name(scale_weight_name) else ""

  keep <- cv_alt$standardizer_upper == std_u &
    cv_alt$detector_upper == det_u &
    cv_alt$T_label == T_label &
    cv_alt$q == q &
    abs(cv_alt$alpha - alpha) < 1e-12

  if (det_u %in% c("FULLCUSUM", "PAGECUSUM", "WEIGHTEDCUSUM")) {
    keep <- keep & abs(cv_alt$gamma - as.numeric(gamma)) < 1e-12
  }
  if (det_u == "MOSUM") {
    keep <- keep & abs(cv_alt$bandwidth_h - as.numeric(bandwidth_h)) < 1e-12
  }
  if (det_u == "WEIGHTEDCUSUM") {
    keep <- keep & cv_alt$omega_norm == omega_norm
  }
  if (det_u == "MULTISCALEMOSUM") {
    keep <- keep & cv_alt$hset_name == hset_name & cv_alt$scale_weight_norm == scale_norm
  }

  sub <- cv_alt[keep, , drop = FALSE]
  if (nrow(sub) != 1L) {
    stop(sprintf("Alternative-detector critical value not uniquely found for (%s,%s,T=%s,q=%s,alpha=%s). Matches=%d",
                 standardizer, detector, T_label, q, alpha, nrow(sub)))
  }
  as.numeric(sub$critical_value[1L])
}
