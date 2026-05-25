# ==============================================================
# dependencies.R -- package checks for the simulation code
# ==============================================================

ensure_packages <- function(pkgs, install_if_missing = FALSE) {
  pkgs <- unique(pkgs[nzchar(pkgs)])
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

simulation_required_packages <- function() c("fdapace", "fda")
plot_required_packages <- function() c("ggplot2")

ensure_simulation_packages <- function(install_if_missing = FALSE) {
  ensure_packages(simulation_required_packages(), install_if_missing = install_if_missing)
}

ensure_plot_packages <- function(install_if_missing = FALSE) {
  ensure_packages(unique(c(simulation_required_packages(), plot_required_packages())), install_if_missing = install_if_missing)
}
