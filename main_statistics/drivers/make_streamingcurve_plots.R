# ==============================================================
# make_streamingcurve_plots.R
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
source(file.path(ROOT, "R", "dependencies.R"), local = FALSE)
source(file.path(ROOT, "R", "utils.R"), local = FALSE)
ensure_plot_packages(install_if_missing = FALSE)

summary_dir <- file.path(ROOT, "outputs", "summary")
plot_dir <- file.path(ROOT, "outputs", "plots")
ensure_dir(plot_dir)

size_path <- file.path(summary_dir, "size_summary.csv")
if (!file.exists(size_path)) stop("Run summarize_streamingcurve_simulations.R first.")
size_df <- read.csv(size_path, stringsAsFactors = FALSE)
if ("message" %in% names(size_df)) stop("size_summary.csv does not contain summary data. Inspect outputs/summary diagnostics.")

power_path <- file.path(summary_dir, "power_summary.csv")
power_df <- if (file.exists(power_path)) read.csv(power_path, stringsAsFactors = FALSE) else data.frame()

library(ggplot2)

if (!"method_label" %in% names(size_df)) size_df$method_label <- size_df$method_id
if (!"method_group" %in% names(size_df)) size_df$method_group <- "All"
size_df$T <- factor(as.character(size_df$T), levels = c("1", "2", "5", "10", "Inf"))

for (grp in unique(size_df$method_group)) {
  for (a in sort(unique(size_df$alpha))) {
    sub <- subset(size_df, method_group == grp & abs(alpha - a) < 1e-12)
    if (nrow(sub) == 0L) next
    p <- ggplot(sub, aes(x = T, y = size, group = method_label, linetype = method_label)) +
      geom_line() + geom_point() +
      facet_grid(dgp_type ~ m, scales = "free_y") +
      labs(x = "T", y = "Empirical size", title = paste0("Null size: ", grp, ", alpha=", a), linetype = "Method") +
      theme_minimal(base_size = 11)
    ggplot2::ggsave(filename = file.path(plot_dir, paste0("size_", sanitize_tag(grp), "_alpha_", sanitize_tag(a), ".png")),
           plot = p, width = 12, height = 8, dpi = 300)
  }
}

if (nrow(power_df) > 0L && !"message" %in% names(power_df)) {
  if (!"method_label" %in% names(power_df)) power_df$method_label <- power_df$method_id
  if (!"method_group" %in% names(power_df)) power_df$method_group <- "All"
  power_df$T <- factor(as.character(power_df$T), levels = c("1", "2", "5", "10", "Inf"))
  power_df$scenario <- factor(power_df$scenario, levels = unique(power_df$scenario))

  for (sc in unique(power_df$scenario)) {
    for (grp in unique(power_df$method_group)) {
      for (a in sort(unique(power_df$alpha))) {
        sub <- subset(power_df, scenario == sc & method_group == grp & abs(alpha - a) < 1e-12)
        if (nrow(sub) == 0L) next

        p1 <- ggplot(sub, aes(x = delta, y = power, group = method_label, linetype = method_label)) +
          geom_line() + geom_point(size = 0.8) +
          facet_grid(dgp_type + s_star ~ T, scales = "free_y") +
          labs(x = expression(Delta), y = "Power", title = paste0("Power: ", sc, ", ", grp, ", alpha=", a), linetype = "Method") +
          theme_minimal(base_size = 11)
        ggplot2::ggsave(filename = file.path(plot_dir, paste0("power_", sanitize_tag(sc), "_", sanitize_tag(grp), "_alpha_", sanitize_tag(a), ".png")),
               plot = p1, width = 14, height = 10, dpi = 300)

        if ("size_adjusted_power" %in% names(sub)) {
          p2 <- ggplot(sub, aes(x = delta, y = size_adjusted_power, group = method_label, linetype = method_label)) +
            geom_line() + geom_point(size = 0.8) +
            facet_grid(dgp_type + s_star ~ T, scales = "free_y") +
            labs(x = expression(Delta), y = "Size-adjusted power", title = paste0("Size-adjusted power: ", sc, ", ", grp, ", alpha=", a), linetype = "Method") +
            theme_minimal(base_size = 11)
          ggplot2::ggsave(filename = file.path(plot_dir, paste0("size_adjusted_power_", sanitize_tag(sc), "_", sanitize_tag(grp), "_alpha_", sanitize_tag(a), ".png")),
                 plot = p2, width = 14, height = 10, dpi = 300)
        }
      }
    }
  }
}

message("Plots written to: ", normalize_path2(plot_dir, mustWork = FALSE))
