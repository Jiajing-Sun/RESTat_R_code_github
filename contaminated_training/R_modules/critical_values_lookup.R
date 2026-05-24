# ==============================================================
# critical_values_lookup.R -- load and query the latest critical values
# ============================================================== 

load_main_critical_values <- function(root) {
  base_path <- file.path(root, "critical_values", "critical_values_all.csv")
  weight_path <- file.path(root, "critical_values", "critical_values_all_weights.csv")
  if (!file.exists(base_path)) stop("Missing main critical value file: ", base_path)
  if (!file.exists(weight_path)) stop("Missing weighted critical value file: ", weight_path)
  base <- read.csv(base_path, stringsAsFactors = FALSE)
  weights <- read.csv(weight_path, stringsAsFactors = FALSE)
  base$T_label <- canonical_T_vec(base$T)
  weights$T_label <- canonical_T_vec(weights$T)
  base$stat_upper <- toupper(base$stat)
  base$type_upper <- toupper(base$type)
  weights$stat_upper <- toupper(weights$stat)
  weights$type_upper <- toupper(weights$type)
  weights$weight_norm <- vapply(weights$weight_name, normalize_weight_name, character(1L))
  list(base = base, weights = weights)
}

load_alt_critical_values <- function(root) {
  path <- file.path(root, "critical_values", "critical_values_alt_detectors.csv")
  if (!file.exists(path)) stop("Missing alternative-detector critical value file: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE)
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
