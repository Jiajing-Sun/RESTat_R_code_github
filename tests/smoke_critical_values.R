args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
root <- if (length(file_arg) > 0L) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L]))
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

out_dir <- file.path(tempdir(), paste0("streaming_curve_cv_smoke_", Sys.getpid()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cv_dir <- file.path(root, "critical_values", "main")
owd <- setwd(cv_dir)
on.exit(setwd(owd), add = TRUE)
cmd <- "run_generate_streamingcurve_cv.R"
args <- c(
  cmd,
  "--q_grid=1",
  "--T_grid=1",
  "--nrep_finite=5",
  "--grid_monitor=20",
  "--ncores=1",
  "--include_openend=false",
  paste0("--output_dir=", out_dir),
  "--backup_existing_outputs=false",
  "--verbose=false"
)

status <- system2("Rscript", args = args)
if (!identical(status, 0L)) {
  quit(status = status)
}

expected <- file.path(out_dir, "critical_values_all.csv")
if (!file.exists(expected)) {
  message("Smoke test did not create expected file: ", expected)
  quit(status = 1L)
}

message("Critical-value smoke test wrote: ", expected)
