# ==============================================================
# make_contaminated_training_plots.R
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
alpha_target <- as.numeric(args$alpha %||% 0.05)
output_tier <- sanitize_tag(args$output_tier %||% "final_results")
output_root <- file.path(ROOT, "outputs", output_tier)
summary_dir <- file.path(output_root, "summary")
plot_dir <- file.path(output_root, "plots_contam_training", scenario)
ensure_dir(plot_dir)
summary_file <- file.path(summary_dir, paste0("contam_training_summary_", scenario, ".csv"))
if (!file.exists(summary_file)) stop("Run summarize_contaminated_training_streamingcurve.R first: missing ", summary_file)

sumdf <- read.csv(summary_file, stringsAsFactors = FALSE)
if (nrow(sumdf) == 0L) stop("Summary file is empty: ", summary_file)
for (nm in c("alpha", "train_b", "break_frac", "post_break_detect_rate", "pre_break_stop_rate", "conditional_delay")) {
  sumdf[[nm]] <- suppressWarnings(as.numeric(sumdf[[nm]]))
}
sub <- subset(sumdf, abs(alpha - alpha_target) < 1e-12)
if (nrow(sub) == 0L) stop("No rows for alpha=", alpha_target)

plot_metric <- function(df, metric, metric_label, outfile_prefix) {
  method_ids <- unique(df$method_id)
  dgp_levels <- sort(unique(df$dgp_type))
  contam_levels <- c("late_break", "drift")
  break_levels <- sort(unique(df$break_frac))
  for (mid in method_ids) {
    d0 <- subset(df, method_id == mid)
    if (nrow(d0) == 0L) next
    label <- unique(na.omit(d0$method_label))[1L]
    if (is.na(label) || !nzchar(label)) label <- mid
    for (bf in break_levels) {
      d1 <- subset(d0, abs(break_frac - bf) < 1e-10)
      if (nrow(d1) == 0L) next
      out_png <- file.path(plot_dir, sprintf("%s_%s_break%s.png", outfile_prefix, sanitize_tag(mid), gsub("\\.", "", format(bf, nsmall = 2))))
      png(out_png, width = 1600, height = 900, res = 140)
      par(mfrow = c(length(contam_levels), length(dgp_levels)), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
      for (ct in contam_levels) {
        for (dg in dgp_levels) {
          dd <- subset(d1, train_contam == ct & dgp_type == dg)
          dd <- dd[order(dd$train_b), , drop = FALSE]
          if (nrow(dd) == 0L) {
            plot.new(); title(main = paste(ct, dg)); next
          }
          x <- dd$train_b; y <- dd[[metric]]
          plot(x, y, type = "o", pch = 19, lwd = 2,
               ylim = c(0, if (metric == "conditional_delay") max(y, na.rm = TRUE) * 1.1 else 1),
               xlab = expression(b[train]), ylab = metric_label,
               main = sprintf("%s | %s", ct, dg))
          grid()
        }
      }
      mtext(sprintf("%s | break fraction = %.2f | alpha = %.2f", label, bf, alpha_target), outer = TRUE, cex = 1.1, font = 2)
      dev.off()
    }
  }
}

plot_metric(sub, "post_break_detect_rate", "Post-break detection", "contam_postbreak")
plot_metric(sub, "pre_break_stop_rate", "Pre-break stop rate", "contam_prebreak")
plot_metric(sub, "conditional_delay", "Conditional delay", "contam_delay")

message("Wrote contaminated-training plots to: ", normalize_path2(plot_dir, mustWork = FALSE))
