# ==============================================================
# critical_values.R -- user-facing loader / lookup helpers
# ============================================================== 

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

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p1)) return(normalizePath(p1, winslash = "/", mustWork = FALSE))
  }

  NULL
}

source_project_paths <- function() {
  p <- bootstrap_script_path()
  candidates <- unique(c(
    if (!is.null(p)) file.path(dirname(p), "project_paths.R") else NULL,
    file.path(getwd(), "R", "project_paths.R"),
    file.path(getwd(), "project_paths.R")
  ))

  for (x in candidates) {
    if (!is.null(x) && file.exists(x)) {
      source(x, local = FALSE)
      return(invisible(x))
    }
  }

  stop("Could not locate project_paths.R while sourcing critical_values.R.")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
source_project_paths()
ROOT <- resolve_project_root(default_start = dirname(bootstrap_script_path() %||% getwd()))

source(file.path(ROOT, "R", "utils.R"), local = FALSE)
source(file.path(ROOT, "R", "weights.R"), local = FALSE)
source(file.path(ROOT, "R", "critical_values_io.R"), local = FALSE)
