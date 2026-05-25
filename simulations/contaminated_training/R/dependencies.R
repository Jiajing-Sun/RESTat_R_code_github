# ==============================================================
# dependencies.R -- minimal package checks for the Windows-safe
# streaming-curve simulation code
# ============================================================== 

ensure_packages <- function(pkgs, install_if_missing = FALSE) {
  pkgs <- unique(pkgs[nzchar(pkgs)])
  if (length(pkgs) == 0L) return(invisible(character(0)))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L && install_if_missing) {
    install.packages(missing, dependencies = TRUE)
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(still_missing) > 0L) {
    stop(
      "Missing required R packages: ", paste(still_missing, collapse = ", "),
      ". Install them first, or rerun with --install_missing=TRUE if appropriate."
    )
  }
  invisible(pkgs)
}

# The simulation engine below is intentionally written in base R + recommended
# packages only, so that it is less fragile on Windows PSOCK clusters.
simulation_required_packages <- function() character(0)
plot_required_packages <- function() character(0)

ensure_simulation_packages <- function(install_if_missing = FALSE) {
  invisible(simulation_required_packages())
}

ensure_plot_packages <- function(install_if_missing = FALSE) {
  invisible(plot_required_packages())
}
