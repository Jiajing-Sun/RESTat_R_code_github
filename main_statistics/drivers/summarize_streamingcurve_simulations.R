# ==============================================================
# summarize_streamingcurve_simulations.R
# ============================================================== 

bootstrap_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p1)) return(normalizePath(p1, winslash = "/", mustWork = FALSE))
    p2 <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  NULL
}
source_project_paths_bootstrap <- function() {
  p <- bootstrap_script_path()
  candidates <- unique(c(
    if (!is.null(p)) file.path(dirname(p), "R", "project_paths.R") else NULL,
    file.path(getwd(), "R", "project_paths.R"),
    file.path(dirname(getwd()), "R", "project_paths.R")
  ))
  for (x in candidates) {
    if (!is.null(x) && file.exists(x)) {
      source(x, local = FALSE)
      return(invisible(x))
    }
  }
  stop("Could not locate R/project_paths.R during bootstrap.")
}
source_project_paths_bootstrap()
ROOT <- resolve_project_root(default_start = dirname(bootstrap_script_path() %||% getwd()))
source(file.path(ROOT, "R", "project_paths.R"), local = FALSE)
source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "method_catalog.R"), local = FALSE)

null_dir <- file.path(ROOT, "outputs", "null_raw")
power_root <- file.path(ROOT, "outputs", "power_raw")
summary_dir <- file.path(ROOT, "outputs", "summary")
ensure_dir(summary_dir)

expected_cols <- c(
  "rep", "dgp_type", "scenario", "m", "T", "delta", "s_star", "q_used", "fpca_mode", "fve_threshold",
  "family", "standardizer", "detector", "type", "gamma", "weight_name", "bandwidth_h",
  "omega_name", "hset_name", "scale_weight_name",
  "method_group", "method_label", "method_id",
  "alpha", "statistic", "critical_value", "reject", "first_rejection",
  "error_flag", "error_message", "source_file"
)

ensure_col <- function(df, nm, value = NA) {
  if (!nm %in% names(df)) df[[nm]] <- value
  df
}

normalize_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  y <- trimws(toupper(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %in% c("FALSE", "F", "0", "NO", "N", "")] <- FALSE
  out
}

canonicalize_T_col <- function(x) {
  out <- rep(NA_character_, length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi)) out[i] <- NA_character_ else out[i] <- tryCatch(canonical_T_scalar(xi), error = function(e) as.character(xi))
  }
  out
}

infer_method_label <- function(df) {
  out <- rep(NA_character_, nrow(df))
  has_main <- !is.na(df$family) & df$family == "Main" & !is.na(df$standardizer) & !is.na(df$type)
  if (any(has_main)) {
    out[has_main] <- mapply(method_label_main,
                            df$standardizer[has_main], df$type[has_main], df$gamma[has_main], df$weight_name[has_main],
                            USE.NAMES = FALSE)
  }
  has_alt <- !is.na(df$detector) & df$detector %in% c("FullCUSUM", "PageCUSUM", "MOSUM", "WeightedCUSUM", "MultiscaleMOSUM")
  if (any(has_alt)) {
    out[has_alt] <- mapply(method_label_alt,
                           df$standardizer[has_alt], df$detector[has_alt], df$gamma[has_alt], df$bandwidth_h[has_alt],
                           df$omega_name[has_alt], df$hset_name[has_alt], df$scale_weight_name[has_alt],
                           USE.NAMES = FALSE)
  }
  out
}

infer_method_group <- function(df) {
  out <- rep(NA_character_, nrow(df))
  ok <- !is.na(df$family) | !is.na(df$type) | !is.na(df$detector)
  out[ok] <- mapply(method_group_from_row,
                    ifelse(is.na(df$family[ok]), "", df$family[ok]),
                    ifelse(is.na(df$type[ok]), "", df$type[ok]),
                    ifelse(is.na(df$detector[ok]), "", df$detector[ok]),
                    USE.NAMES = FALSE)
  out
}

infer_method_id <- function(df) {
  base <- ifelse(!is.na(df$method_label) & nzchar(df$method_label), df$method_label,
                 paste(df$family, df$standardizer, df$detector, df$type,
                       ifelse(is.na(df$gamma), "", paste0("g", df$gamma)),
                       ifelse(is.na(df$weight_name), "", df$weight_name),
                       ifelse(is.na(df$bandwidth_h), "", paste0("h", df$bandwidth_h)),
                       ifelse(is.na(df$omega_name), "", df$omega_name),
                       ifelse(is.na(df$hset_name), "", df$hset_name),
                       ifelse(is.na(df$scale_weight_name), "", df$scale_weight_name)))
  sanitize_tag(tolower(gsub("[^A-Za-z0-9]+", "_", base)))
}

normalize_raw_df <- function(df, source_file) {
  if (nrow(df) == 0L) {
    out <- as.data.frame(setNames(replicate(length(expected_cols), logical(0), simplify = FALSE), expected_cols))
    return(out)
  }

  renames <- c(stat = "statistic", cv = "critical_value", reject_flag = "reject",
               first_stop = "first_rejection", std = "standardizer")
  for (old_nm in names(renames)) {
    new_nm <- renames[[old_nm]]
    if (old_nm %in% names(df) && !new_nm %in% names(df)) names(df)[names(df) == old_nm] <- new_nm
  }

  for (nm in expected_cols) df <- ensure_col(df, nm, NA)
  df$source_file <- source_file
  df$error_flag <- normalize_bool(df$error_flag)
  if (all(is.na(df$error_flag))) {
    err_msg_chr <- ifelse(is.na(df$error_message), "", as.character(df$error_message))
    df$error_flag <- nzchar(err_msg_chr) & is.na(df$statistic) & is.na(df$critical_value) & is.na(df$reject)
  }

  df$T <- canonicalize_T_col(df$T)
  num_cols <- c("rep","m","delta","s_star","q_used","fve_threshold","gamma","bandwidth_h",
                "alpha","statistic","critical_value","first_rejection")
  for (nm in num_cols) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  df$reject <- normalize_bool(df$reject)

  chr_cols <- c("dgp_type","scenario","family","standardizer","detector","type","weight_name","omega_name",
                "hset_name","scale_weight_name","method_group","method_label","method_id","error_message","fpca_mode")
  for (nm in chr_cols) df[[nm]] <- ifelse(is.na(df[[nm]]), NA_character_, as.character(df[[nm]]))

  need_label <- is.na(df$method_label) | !nzchar(df$method_label)
  if (any(need_label)) df$method_label[need_label] <- infer_method_label(df)[need_label]

  need_group <- is.na(df$method_group) | !nzchar(df$method_group)
  if (any(need_group)) df$method_group[need_group] <- infer_method_group(df)[need_group]

  need_id <- is.na(df$method_id) | !nzchar(df$method_id)
  if (any(need_id)) df$method_id[need_id] <- infer_method_id(df)[need_id]

  for (nm in c("weight_name","omega_name","hset_name","scale_weight_name")) {
    df[[nm]][!nzchar(ifelse(is.na(df[[nm]]), "", df[[nm]]))] <- NA_character_
  }

  df[, expected_cols, drop = FALSE]
}

safe_read_one_csv <- function(f) {
  tryCatch({
    df <- read.csv(f, stringsAsFactors = FALSE)
    normalize_raw_df(df, basename(f))
  }, error = function(e) {
    data.frame(
      rep = NA_integer_, dgp_type = NA_character_, scenario = NA_character_, m = NA_real_, T = NA_character_,
      delta = NA_real_, s_star = NA_real_, q_used = NA_real_, fpca_mode = NA_character_, fve_threshold = NA_real_,
      family = NA_character_, standardizer = NA_character_, detector = NA_character_, type = NA_character_,
      gamma = NA_real_, weight_name = NA_character_, bandwidth_h = NA_real_, omega_name = NA_character_,
      hset_name = NA_character_, scale_weight_name = NA_character_,
      method_group = "ReadError", method_label = "ReadError", method_id = NA_character_,
      alpha = NA_real_, statistic = NA_real_, critical_value = NA_real_, reject = NA, first_rejection = NA_real_,
      error_flag = TRUE, error_message = paste0("Failed to read CSV: ", conditionMessage(e)), source_file = basename(f),
      stringsAsFactors = FALSE
    )
  })
}

read_many_csv <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0L) return(data.frame())
  out <- lapply(files, safe_read_one_csv)
  do.call(rbind, out)
}

fill_group_nas <- function(df) {
  chr_cols <- names(df)[vapply(df, is.character, logical(1L))]
  for (nm in chr_cols) df[[nm]][is.na(df[[nm]])] <- ""
  num_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  for (nm in num_cols) df[[nm]][is.na(df[[nm]])] <- -999999
  df
}
restore_group_nas <- function(df) {
  chr_cols <- names(df)[vapply(df, is.character, logical(1L))]
  for (nm in chr_cols) df[[nm]][df[[nm]] == ""] <- NA_character_
  num_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  for (nm in num_cols) df[[nm]][df[[nm]] == -999999] <- NA_real_
  df
}

write_diag_table <- function(df, path) {
  write.csv(df, path, row.names = FALSE)
}

write_raw_diagnostics <- function(df, output_prefix) {
  if (nrow(df) == 0L) {
    write_diag_table(data.frame(message = "No raw CSV files found."), file.path(summary_dir, paste0(output_prefix, "_raw_diagnostics.csv")))
    return(invisible(NULL))
  }
  diag_df <- aggregate(rep(1L, nrow(df)) ~ source_file + error_flag, data = fill_group_nas(df), FUN = length)
  names(diag_df)[names(diag_df) == "rep(1L, nrow(df))"] <- "n_rows"
  diag_df <- restore_group_nas(diag_df)
  write_diag_table(diag_df, file.path(summary_dir, paste0(output_prefix, "_raw_diagnostics.csv")))

  err_df <- subset(df, isTRUE(error_flag) | (!is.na(error_flag) & error_flag))
  if (nrow(err_df) > 0L) {
    err_df$error_message <- ifelse(is.na(err_df$error_message) | !nzchar(err_df$error_message), "<missing>", err_df$error_message)
    err_sum <- aggregate(rep(1L, nrow(err_df)) ~ error_message, data = err_df, FUN = length)
    names(err_sum)[2] <- "n_rows"
    err_sum <- err_sum[order(-err_sum$n_rows), , drop = FALSE]
    write_diag_table(err_sum, file.path(summary_dir, paste0(output_prefix, "_raw_error_messages.csv")))
  }
  invisible(NULL)
}

build_size_summary <- function(usable_null) {
  usable_null$reject_num <- as.numeric(usable_null$reject)
  grp <- aggregate(cbind(reject_num, rep) ~ dgp_type + m + T + alpha + method_group + method_id + method_label + family + standardizer + detector + type + gamma + weight_name + bandwidth_h + omega_name + hset_name + scale_weight_name,
                   data = fill_group_nas(usable_null),
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE), n = sum(!is.na(x))))
  out <- grp[, setdiff(names(grp), c("reject_num", "rep")), drop = FALSE]
  out$size <- grp$reject_num[, "mean"]
  out$n_valid <- grp$reject_num[, "n"]
  restore_group_nas(out)
}

build_empirical_null_cv <- function(usable_null) {
  key <- unique(usable_null[, c("dgp_type", "m", "T", "alpha", "method_id")])
  key$empirical_critical_value <- NA_real_
  for (i in seq_len(nrow(key))) {
    sub <- usable_null$statistic[
      usable_null$dgp_type == key$dgp_type[i] &
      usable_null$m == key$m[i] &
      usable_null$T == key$T[i] &
      abs(usable_null$alpha - key$alpha[i]) < 1e-12 &
      usable_null$method_id == key$method_id[i]
    ]
    key$empirical_critical_value[i] <- as.numeric(stats::quantile(sub, probs = 1 - key$alpha[i], na.rm = TRUE, names = FALSE, type = 7))
  }
  key
}

build_power_summary <- function(usable_power, emp_null_cv) {
  usable_power$reject_num <- as.numeric(usable_power$reject)
  power_df2 <- merge(usable_power, emp_null_cv, by = c("dgp_type", "m", "T", "alpha", "method_id"), all.x = TRUE, sort = FALSE)
  power_df2$size_adjusted_reject <- with(power_df2, as.numeric(!is.na(empirical_critical_value) & statistic > empirical_critical_value))
  grp <- aggregate(cbind(reject_num, size_adjusted_reject) ~ scenario + dgp_type + m + T + s_star + delta + alpha + method_group + method_id + method_label + family + standardizer + detector + type + gamma + weight_name + bandwidth_h + omega_name + hset_name + scale_weight_name,
                   data = fill_group_nas(power_df2),
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE), n = sum(!is.na(x))))
  out <- grp[, setdiff(names(grp), c("reject_num", "size_adjusted_reject")), drop = FALSE]
  out$power <- grp$reject_num[, "mean"]
  out$size_adjusted_power <- grp$size_adjusted_reject[, "mean"]
  out$n_valid <- grp$reject_num[, "n"]
  restore_group_nas(out)
}

null_df <- read_many_csv(null_dir)
write_raw_diagnostics(null_df, "null")
if (nrow(null_df) == 0L) stop("No null raw simulation files found in ", null_dir)

usable_null <- subset(null_df,
                      (is.na(error_flag) | !error_flag) &
                      !is.na(method_id) & nzchar(method_id) &
                      !is.na(alpha) & !is.na(statistic) & !is.na(reject))

if (nrow(usable_null) == 0L) {
  stop(paste0(
    "Null raw CSVs were found, but no usable simulation rows remained after filtering. ",
    "Most likely the raw files contain only error rows. ",
    "Please inspect outputs/summary/null_raw_diagnostics.csv and outputs/summary/null_raw_error_messages.csv."
  ))
}

size_summary <- build_size_summary(usable_null)
write.csv(size_summary, file.path(summary_dir, "size_summary.csv"), row.names = FALSE)

emp_null_cv <- build_empirical_null_cv(usable_null)
write.csv(emp_null_cv, file.path(summary_dir, "empirical_null_critical_values.csv"), row.names = FALSE)

power_dirs <- list.dirs(power_root, recursive = FALSE, full.names = TRUE)
if (length(power_dirs) > 0L) {
  power_df <- do.call(rbind, lapply(power_dirs, read_many_csv))
  write_raw_diagnostics(power_df, "power")
  usable_power <- subset(power_df,
                         (is.na(error_flag) | !error_flag) &
                         !is.na(method_id) & nzchar(method_id) &
                         !is.na(alpha) & !is.na(statistic) & !is.na(reject))
  if (nrow(usable_power) > 0L) {
    power_summary <- build_power_summary(usable_power, emp_null_cv)
    write.csv(power_summary, file.path(summary_dir, "power_summary.csv"), row.names = FALSE)
  } else {
    write.csv(data.frame(message = "Power raw CSVs were found, but no usable rows remained after filtering."),
              file.path(summary_dir, "power_summary.csv"), row.names = FALSE)
  }
}

message("Summary files written to: ", normalize_path2(summary_dir, mustWork = FALSE))
