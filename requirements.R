required_packages <- c(
  "data.table",
  "fda",
  "fdapace",
  "ggplot2",
  "gtable"
)

optional_packages <- c("rstudioapi")

args <- commandArgs(trailingOnly = TRUE)
install_missing <- "--install" %in% args

check_packages <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
}

missing_required <- check_packages(required_packages)
missing_optional <- check_packages(optional_packages)

if (length(missing_required) > 0L && install_missing) {
  install.packages(missing_required, repos = "https://cloud.r-project.org")
  missing_required <- check_packages(required_packages)
}

if (length(missing_required) > 0L) {
  message("Missing required R package(s): ", paste(missing_required, collapse = ", "))
  message("Install them with: Rscript requirements.R --install")
  quit(status = 1L)
}

message("Required R packages are available: ", paste(required_packages, collapse = ", "))
if (length(missing_optional) > 0L) {
  message("Optional R package(s) not installed: ", paste(missing_optional, collapse = ", "))
}
