# ==============================================================
# critical_values_io.R -- clean load/save/lookup helpers
# ============================================================== 

canonicalize_base_table <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(empty_base_table())
  out <- df

  need <- c("stat", "type", "T", "gamma", "q", "alpha", "critical_value", "source_file")
  miss <- setdiff(need, names(out))
  for (nm in miss) {
    if (nm %in% c("gamma", "alpha", "critical_value")) out[[nm]] <- NA_real_
    else if (nm == "q") out[[nm]] <- NA_integer_
    else out[[nm]] <- ""
  }

  out$stat <- as.character(out$stat)
  out$type <- as.character(out$type)
  out$T <- canonical_T_vec(out$T)
  out$gamma <- suppressWarnings(as.numeric(out$gamma))
  out$q <- as.integer(out$q)
  out$alpha <- as.numeric(out$alpha)
  out$critical_value <- as.numeric(out$critical_value)
  out$source_file <- as.character(out$source_file)

  out <- out[, need, drop = FALSE]
  out
}

canonicalize_weight_table <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(empty_weight_table())
  out <- df

  need <- c("stat", "type", "T", "gamma", "q", "alpha", "critical_value", "source_file", "weight_name")
  miss <- setdiff(need, names(out))
  for (nm in miss) {
    if (nm %in% c("gamma", "alpha", "critical_value")) out[[nm]] <- NA_real_
    else if (nm == "q") out[[nm]] <- NA_integer_
    else out[[nm]] <- ""
  }

  out$stat <- as.character(out$stat)
  out$type <- as.character(out$type)
  out$T <- canonical_T_vec(out$T)
  out$gamma <- suppressWarnings(as.numeric(out$gamma))
  out$q <- as.integer(out$q)
  out$alpha <- as.numeric(out$alpha)
  out$critical_value <- as.numeric(out$critical_value)
  out$source_file <- as.character(out$source_file)
  out$weight_name <- vapply(out$weight_name, normalize_weight_name, character(1L))

  out <- out[, need, drop = FALSE]
  out
}

base_key <- function(df) {
  paste0(
    toupper(df$stat), "|",
    toupper(df$type), "|",
    df$T, "|",
    ifelse(is.na(df$gamma), "NA", format(df$gamma, digits = 15)), "|",
    df$q, "|",
    format(df$alpha, digits = 15)
  )
}

weight_key <- function(df) {
  paste0(
    toupper(df$stat), "|",
    toupper(df$type), "|",
    df$T, "|",
    ifelse(is.na(df$gamma), "NA", format(df$gamma, digits = 15)), "|",
    df$q, "|",
    format(df$alpha, digits = 15), "|",
    toupper(df$weight_name)
  )
}

deduplicate_base <- function(df) {
  df <- canonicalize_base_table(df)
  if (nrow(df) == 0L) return(df)
  df[!duplicated(base_key(df), fromLast = TRUE), , drop = FALSE]
}

deduplicate_weights <- function(df) {
  df <- canonicalize_weight_table(df)
  if (nrow(df) == 0L) return(df)
  df[!duplicated(weight_key(df), fromLast = TRUE), , drop = FALSE]
}

append_cv_rows <- function(cv, base_rows_new = NULL, weight_rows_new = NULL) {
  base_all <- deduplicate_base(rbind(canonicalize_base_table(cv$base), canonicalize_base_table(base_rows_new)))
  w_all <- deduplicate_weights(rbind(canonicalize_weight_table(cv$weights), canonicalize_weight_table(weight_rows_new)))
  list(base = base_all, weights = w_all)
}

load_critical_values <- function(base_path = NULL,
                                 weights_path = NULL,
                                 project_root = NULL,
                                 output_dir = "outputs") {
  root <- project_root %||% getwd()
  base <- empty_base_table()
  weights <- empty_weight_table()

  if (is.null(base_path)) base_path <- file.path(root, output_dir, "critical_values_all.csv")
  if (is.null(weights_path)) weights_path <- file.path(root, output_dir, "critical_values_all_weights.csv")

  if (!is.null(base_path) && file.exists(base_path)) {
    base <- deduplicate_base(read.csv(base_path, stringsAsFactors = FALSE))
  }
  if (!is.null(weights_path) && file.exists(weights_path)) {
    weights <- deduplicate_weights(read.csv(weights_path, stringsAsFactors = FALSE))
  }

  list(base = base, weights = weights)
}

save_critical_values <- function(cv,
                                 out_dir,
                                 base_name = "critical_values_all.csv",
                                 weights_name = "critical_values_all_weights.csv") {
  ensure_dir(out_dir)
  base_path <- file.path(out_dir, base_name)
  weights_path <- file.path(out_dir, weights_name)

  utils::write.csv(deduplicate_base(cv$base), base_path, row.names = FALSE)
  utils::write.csv(deduplicate_weights(cv$weights), weights_path, row.names = FALSE)

  invisible(list(base_path = base_path, weights_path = weights_path))
}

prepare_fresh_output_files <- function(out_dir,
                                       base_name = "critical_values_all.csv",
                                       weights_name = "critical_values_all_weights.csv",
                                       backup_existing = TRUE) {
  ensure_dir(out_dir)
  base_path <- file.path(out_dir, base_name)
  weights_path <- file.path(out_dir, weights_name)

  if (backup_existing) {
    backup_file_if_exists(base_path)
    backup_file_if_exists(weights_path)
  }

  if (file.exists(base_path)) file.remove(base_path)
  if (file.exists(weights_path)) file.remove(weights_path)

  invisible(list(base_path = base_path, weights_path = weights_path))
}

get_critical_value <- function(cv,
                               stat = c("SSMS", "RSMS", "HAC"),
                               type = c("KS", "CvM"),
                               T,
                               q,
                               alpha = 0.05,
                               gamma = NA_real_,
                               weight = "U") {
  stat <- toupper(as.character(stat)[1L])
  type <- toupper(as.character(type)[1L])
  Tlab <- canonical_T_scalar(T)
  q <- as.integer(q)
  alpha <- as.numeric(alpha)

  if (type == "KS") {
    df <- canonicalize_base_table(cv$base)
    keep <- (toupper(df$stat) == stat) &
      (toupper(df$type) == "KS") &
      (df$T == Tlab) &
      (df$q == q) &
      (abs(df$alpha - alpha) < 1e-12)

    if (!is.na(gamma)) keep <- keep & row_match_numeric(df$gamma, gamma)
    sub <- df[keep, , drop = FALSE]

    if (nrow(sub) != 1L) {
      stop(sprintf(
        "Critical value not uniquely found for (%s, KS, T=%s, q=%s, alpha=%s, gamma=%s). Matches=%s.",
        stat, Tlab, q, alpha, gamma, nrow(sub)
      ))
    }
    return(sub$critical_value[1L])
  }

  if (type == "CVM") {
    w <- normalize_weight_name(weight)

    if (w == "U") {
      df <- canonicalize_base_table(cv$base)
      keep <- (toupper(df$stat) == stat) &
        (toupper(df$type) == "CVM") &
        (df$T == Tlab) &
        (df$q == q) &
        (abs(df$alpha - alpha) < 1e-12)

      if (!is.na(gamma)) keep <- keep & row_match_numeric(df$gamma, gamma)
      sub <- df[keep, , drop = FALSE]

      if (nrow(sub) != 1L) {
        stop(sprintf(
          "CvM critical value (U) not uniquely found for (%s, T=%s, q=%s, alpha=%s, gamma=%s). Matches=%s.",
          stat, Tlab, q, alpha, gamma, nrow(sub)
        ))
      }
      return(sub$critical_value[1L])
    }

    df <- canonicalize_weight_table(cv$weights)
    keep <- (toupper(df$stat) == stat) &
      (toupper(df$type) == "CVM") &
      (df$T == Tlab) &
      (df$q == q) &
      (abs(df$alpha - alpha) < 1e-12) &
      (toupper(df$weight_name) == toupper(w))

    if (!is.na(gamma)) keep <- keep & row_match_numeric(df$gamma, gamma)
    sub <- df[keep, , drop = FALSE]

    if (nrow(sub) != 1L) {
      stop(sprintf(
        "CvM critical value (%s) not uniquely found for (%s, T=%s, q=%s, alpha=%s, gamma=%s). Matches=%s.",
        w, stat, Tlab, q, alpha, gamma, nrow(sub)
      ))
    }
    return(sub$critical_value[1L])
  }

  stop("Unknown type: ", type, ". Use 'KS' or 'CvM'.")
}
