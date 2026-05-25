decode_rscript_path <- function(x) {
  gsub("~\\+~", " ", x)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(decode_rscript_path(sub("^--file=", "", file_arg[1])))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  named <- grep("^--[A-Za-z0-9_.-]+=", args, value = TRUE)
  for (x in named) {
    key <- sub("^--([^=]+)=.*$", "\\1", x)
    val <- sub("^--[^=]+=", "", x)
    out[[key]] <- val
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(y)
  x
}

run_step <- function(script_dir, script_name, required_files = character()) {
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0L) {
    message("Skipping ", script_name, " because required file(s) are missing:")
    message("  ", paste(missing, collapse = "\n  "))
    return(invisible(FALSE))
  }

  message("Running ", script_name, " ...")
  ok <- tryCatch({
    source(file.path(script_dir, script_name), local = TRUE)
    TRUE
  }, error = function(e) {
    message("Skipping ", script_name, " after error: ", conditionMessage(e))
    FALSE
  })
  invisible(ok)
}

args <- parse_named_args()
script_dir <- get_script_dir()
data_file <- args$data_file %||% Sys.getenv("SPX_DATA_FILE", unset = file.path(script_dir, "data", "SPX.csv"))
Sys.setenv(SPX_DATA_FILE = data_file)

saved_main <- file.path(script_dir, "data", "empirical_main_results.csv")
saved_alt <- file.path(script_dir, "data", "empirical_alt_results.csv")

status <- c(
  spx_summary = run_step(script_dir, "prepare_spx_summary_figures.R", data_file),
  curve_regimes = run_step(script_dir, "build_intraday_curve_regime_figure.R", data_file),
  surface_detection = run_step(script_dir, "build_intraday_surface_detection_figure.R", data_file),
  saved_result_tables = run_step(
    script_dir,
    "rebuild_empirical_outputs_from_saved_results.R",
    c(saved_main, saved_alt)
  ),
  decision_window = run_step(script_dir, "compute_decision_window_diagnostics.R", data_file),
  q_robustness = run_step(script_dir, "compute_q_robustness_diagnostics.R", data_file),
  training_window = run_step(script_dir, "compute_training_window_robustness.R", data_file)
)

message("Empirical bundle refresh finished.")
message("Completed steps: ", paste(names(status)[status], collapse = ", "))
if (any(!status)) {
  message("Skipped steps: ", paste(names(status)[!status], collapse = ", "))
}
