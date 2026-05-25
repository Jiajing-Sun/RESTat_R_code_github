# ==============================================================
# critical_values_io_alt.R -- load/save/lookup helpers
# ============================================================== 

canonicalize_alt_table <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(empty_alt_table())
  out <- df

  need <- c("standardizer", "detector", "T", "gamma", "q", "alpha", "critical_value",
            "source_file", "bandwidth_h", "omega_name", "hset_name", "scale_weight_name", "exact_scan")
  miss <- setdiff(need, names(out))
  for (nm in miss) {
    if (nm %in% c("gamma", "alpha", "critical_value", "bandwidth_h")) out[[nm]] <- NA_real_
    else if (nm == "q") out[[nm]] <- NA_integer_
    else if (nm == "exact_scan") out[[nm]] <- FALSE
    else out[[nm]] <- ""
  }

  out$standardizer <- toupper(as.character(out$standardizer))
  out$detector <- as.character(out$detector)
  out$T <- canonical_T_vec(out$T)
  out$gamma <- suppressWarnings(as.numeric(out$gamma))
  out$q <- as.integer(out$q)
  out$alpha <- as.numeric(out$alpha)
  out$critical_value <- as.numeric(out$critical_value)
  out$source_file <- as.character(out$source_file)
  out$bandwidth_h <- suppressWarnings(as.numeric(out$bandwidth_h))
  omega_chr <- as.character(out$omega_name)
  out$omega_name <- ""
  keep_omega <- !is.na(omega_chr) & nzchar(omega_chr)
  if (any(keep_omega)) out$omega_name[keep_omega] <- vapply(omega_chr[keep_omega], normalize_omega_name, character(1L))
  out$hset_name <- as.character(out$hset_name)
  scale_chr <- as.character(out$scale_weight_name)
  out$scale_weight_name <- ""
  keep_scale <- !is.na(scale_chr) & nzchar(scale_chr)
  if (any(keep_scale)) out$scale_weight_name[keep_scale] <- vapply(scale_chr[keep_scale], normalize_scale_weight_name, character(1L))
  out$exact_scan <- as.logical(out$exact_scan)
  out <- out[, need, drop = FALSE]
  out
}

alt_key <- function(df) {
  paste0(
    toupper(df$standardizer), "|",
    toupper(df$detector), "|",
    df$T, "|",
    ifelse(is.na(df$gamma), "NA", format(df$gamma, digits = 15)), "|",
    df$q, "|",
    format(df$alpha, digits = 15), "|",
    ifelse(is.na(df$bandwidth_h), "NA", format(df$bandwidth_h, digits = 15)), "|",
    toupper(df$omega_name), "|",
    df$hset_name, "|",
    toupper(df$scale_weight_name)
  )
}

deduplicate_alt <- function(df) {
  df <- canonicalize_alt_table(df)
  if (nrow(df) == 0L) return(df)
  df[!duplicated(alt_key(df), fromLast = TRUE), , drop = FALSE]
}

append_alt_rows <- function(cv, rows_new = NULL) {
  deduplicate_alt(rbind(canonicalize_alt_table(cv), canonicalize_alt_table(rows_new)))
}

load_alt_critical_values <- function(path = NULL,
                                     project_root = NULL,
                                     output_dir = "outputs") {
  root <- project_root %||% getwd()
  out <- empty_alt_table()
  if (is.null(path)) path <- file.path(root, output_dir, "critical_values_alt_detectors.csv")
  if (!is.null(path) && file.exists(path)) {
    out <- deduplicate_alt(read.csv(path, stringsAsFactors = FALSE))
  }
  out
}

save_alt_critical_values <- function(cv,
                                     out_dir,
                                     filename = "critical_values_alt_detectors.csv") {
  ensure_dir(out_dir)
  path <- file.path(out_dir, filename)
  utils::write.csv(deduplicate_alt(cv), path, row.names = FALSE)
  invisible(list(path = path))
}

prepare_fresh_alt_output_file <- function(out_dir,
                                          filename = "critical_values_alt_detectors.csv",
                                          backup_existing = TRUE) {
  ensure_dir(out_dir)
  path <- file.path(out_dir, filename)
  if (backup_existing) backup_file_if_exists(path)
  if (file.exists(path)) file.remove(path)
  invisible(list(path = path))
}

get_alt_critical_value <- function(cv,
                                   standardizer = c("HAC", "SSMS", "RSMS"),
                                   detector,
                                   T,
                                   q,
                                   alpha = 0.05,
                                   gamma = NA_real_,
                                   bandwidth_h = NA_real_,
                                   omega_name = "",
                                   hset_name = "",
                                   scale_weight_name = "") {
  df <- canonicalize_alt_table(cv)
  standardizer <- toupper(as.character(standardizer)[1L])
  detector <- as.character(detector)[1L]
  Tlab <- canonical_T_scalar(T)
  q <- as.integer(q)
  alpha <- as.numeric(alpha)
  omega_name <- if (nzchar(omega_name)) normalize_omega_name(omega_name) else ""
  scale_weight_name <- if (nzchar(scale_weight_name)) normalize_scale_weight_name(scale_weight_name) else ""

  keep <- (toupper(df$standardizer) == standardizer) &
    (df$detector == detector) &
    (df$T == Tlab) &
    (df$q == q) &
    (abs(df$alpha - alpha) < 1e-12)

  if (!is.na(gamma)) keep <- keep & row_match_numeric(df$gamma, gamma)
  if (!is.na(bandwidth_h)) keep <- keep & row_match_numeric(df$bandwidth_h, bandwidth_h)
  if (nzchar(omega_name)) keep <- keep & (toupper(df$omega_name) == toupper(omega_name))
  if (nzchar(hset_name)) keep <- keep & (df$hset_name == hset_name)
  if (nzchar(scale_weight_name)) keep <- keep & (toupper(df$scale_weight_name) == toupper(scale_weight_name))

  sub <- df[keep, , drop = FALSE]
  if (nrow(sub) != 1L) {
    stop(sprintf("Critical value not uniquely found for (%s,%s,T=%s,q=%s,alpha=%s). Matches=%s.",
                 standardizer, detector, Tlab, q, alpha, nrow(sub)))
  }
  sub$critical_value[1L]
}
