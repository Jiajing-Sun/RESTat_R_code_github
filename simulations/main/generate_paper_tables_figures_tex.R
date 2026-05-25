suppressPackageStartupMessages({})

bootstrap_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) {
    return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  }
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_path <- bootstrap_script_path()
sim_root <- dirname(script_path)
paper_root <- dirname(sim_root)
main_builder <- file.path(paper_root, "build_latest_simulation_section.R")

if (!file.exists(main_builder)) {
  stop("Could not find build_latest_simulation_section.R at: ", main_builder)
}

Sys.setenv(STREAMINGCURVE_PAPER_ROOT = paper_root)
source(main_builder, local = globalenv(), chdir = FALSE)
