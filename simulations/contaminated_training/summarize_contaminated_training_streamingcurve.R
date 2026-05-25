# ==============================================================
# summarize_contaminated_training_streamingcurve.R
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

args <- parse_named_args()
scenario <- sanitize_tag(args$scenario %||% "level_shift")
output_tier <- sanitize_tag(args$output_tier %||% "final_results")
output_root <- file.path(ROOT, "outputs", output_tier)
raw_root <- file.path(output_root, "raw", scenario)
summary_dir <- file.path(output_root, "summary")
ensure_dir(summary_dir)
include_pattern <- args$include_pattern %||% ""
exclude_pattern <- args$exclude_pattern %||% ""

read_many <- function(root_dir) {
  if (!dir.exists(root_dir)) return(data.frame())
  files_rds <- list.files(root_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  files_csv <- list.files(root_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  files_rds <- files_rds[!grepl("^\\._", basename(files_rds))]
  files_csv <- files_csv[!grepl("^\\._", basename(files_csv))]
  files <- if (length(files_rds) > 0L) files_rds else files_csv
  if (nzchar(include_pattern)) files <- files[grepl(include_pattern, basename(files))]
  if (nzchar(exclude_pattern)) files <- files[!grepl(exclude_pattern, basename(files))]
  if (length(files) == 0L) return(data.frame())
  out <- lapply(files, function(f) {
    ext <- tolower(tools::file_ext(f))
    obj <- if (ext == "rds") safe_read_rds(f) else safe_read_csv(f)
    if (is.null(obj)) return(NULL)
    if (!is.data.frame(obj)) obj <- as.data.frame(obj)
    if (nrow(obj) == 0L) return(NULL)
    obj$source_file <- basename(f)
    obj
  })
  combine_rbind(out)
}

df <- read_many(raw_root)
if (nrow(df) == 0L) stop("No contaminated-training raw outputs were found in: ", raw_root)

for (nm in c("error_flag", "reject")) {
  if (nm %in% names(df)) {
    if (!is.logical(df[[nm]])) {
      y <- trimws(toupper(as.character(df[[nm]])))
      out <- rep(NA, length(y))
      out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
      out[y %in% c("FALSE", "F", "0", "NO", "N", "")] <- FALSE
      df[[nm]] <- out
    }
  } else {
    df[[nm]] <- FALSE
  }
}

# Basic normalization
for (nm in c("m", "T", "delta", "s_star", "q_used", "gamma", "bandwidth_h", "alpha", "statistic", "critical_value", "first_rejection", "train_b", "train_break_frac", "train_drift_start_frac", "t_train_star")) {
  if (nm %in% names(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
}
for (nm in c("dgp_type", "scenario", "family", "standardizer", "detector", "type", "weight_name", "omega_name", "hset_name", "scale_weight_name", "method_group", "method_label", "method_id", "train_contam", "contam_shape")) {
  if (!nm %in% names(df)) df[[nm]] <- NA_character_ else df[[nm]] <- as.character(df[[nm]])
}

usable <- subset(df, !isTRUE(error_flag) & !is.na(method_id) & nzchar(method_id))
if (nrow(usable) == 0L) {
  atomic_write_csv(df, file.path(summary_dir, paste0("contam_training_all_rows_", scenario, ".csv")))
  stop("Contaminated-training raw outputs were found, but no usable rows remained after filtering.")
}

usable$break_frac <- with(usable, ifelse(is.finite(m) & is.finite(T) & is.finite(s_star) & (m * T) > 0, s_star / (m * T), NA_real_))
usable$post_break_detect <- with(usable, !is.na(first_rejection) & !is.na(s_star) & (first_rejection > s_star))
usable$pre_break_stop <- with(usable, !is.na(first_rejection) & !is.na(s_star) & (first_rejection <= s_star))
usable$delay_if_post <- with(usable, ifelse(post_break_detect, first_rejection - s_star, NA_real_))

key_cols <- c("dgp_type", "scenario", "m", "T", "delta", "break_frac", "s_star",
              "train_contam", "train_b", "train_break_frac", "train_drift_start_frac", "contam_shape",
              "method_group", "method_id", "method_label", "family", "standardizer", "detector", "type",
              "gamma", "weight_name", "bandwidth_h", "omega_name", "hset_name", "scale_weight_name", "alpha")

key_df <- usable[, key_cols, drop = FALSE]
for (nm in names(key_df)) {
  x <- key_df[[nm]]
  if (is.numeric(x)) {
    x_chr <- ifelse(is.na(x), "__NA__", format(x, digits = 15, trim = TRUE, scientific = FALSE))
  } else {
    x_chr <- trimws(as.character(x))
    x_chr[is.na(x_chr) | !nzchar(x_chr)] <- "__NA__"
  }
  key_df[[nm]] <- x_chr
}
split_key <- interaction(key_df, drop = TRUE, lex.order = TRUE)
out_list <- lapply(split(usable, split_key), function(dd) {
  base <- dd[1L, key_cols, drop = FALSE]
  base$reject_rate <- mean(dd$reject, na.rm = TRUE)
  base$post_break_detect_rate <- mean(dd$post_break_detect, na.rm = TRUE)
  base$pre_break_stop_rate <- mean(dd$pre_break_stop, na.rm = TRUE)
  base$conditional_delay <- if (all(is.na(dd$delay_if_post))) NA_real_ else mean(dd$delay_if_post, na.rm = TRUE)
  base$n_rep <- nrow(dd)
  base$n_error_rows_ignored <- sum(is.na(dd$reject))
  base
})
summary_df <- combine_rbind(out_list)
summary_df <- summary_df[order(summary_df$method_label, summary_df$dgp_type, summary_df$train_contam, summary_df$train_b, summary_df$break_frac), , drop = FALSE]

atomic_write_csv(summary_df, file.path(summary_dir, paste0("contam_training_summary_", scenario, ".csv")))
atomic_write_csv(usable, file.path(summary_dir, paste0("contam_training_usable_rows_", scenario, ".csv")))

message("Wrote contaminated-training summary to: ", file.path(summary_dir, paste0("contam_training_summary_", scenario, ".csv")))
