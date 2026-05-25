args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
root <- if (length(file_arg) > 0L) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L]))
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

out_file <- file.path(tempdir(), paste0("streaming_curve_alt_cv_smoke_", Sys.getpid(), ".csv"))
cv_dir <- file.path(root, "critical_values", "alternative_detectors")
owd <- setwd(cv_dir)
on.exit(setwd(owd), add = TRUE)

args <- c(
  "generate_alt_detector_critical_values.R",
  "--q_grid=1",
  "--T_grid=1",
  "--finite_nrep=5",
  "--finite_kmax=20",
  "--page_length_grid_size=5",
  "--weighted_length_grid_size=5",
  "--ncores=1",
  "--include_openend=false",
  paste0("--output=", out_file)
)

status <- system2("Rscript", args = args)
if (!identical(status, 0L)) {
  quit(status = status)
}

if (!file.exists(out_file)) {
  message("Alternative-detector smoke test did not create expected file: ", out_file)
  quit(status = 1L)
}

message("Alternative-detector smoke test wrote: ", out_file)
