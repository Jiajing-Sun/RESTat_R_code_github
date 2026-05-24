# ==============================================================
# build_contaminated_training_paper_assets.R
# Build tables, figures, and notes for the targeted
# contaminated-training robustness study.
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
scenario_target <- sanitize_tag(args$scenario %||% "level_shift")
alpha_target <- as.numeric(args$alpha %||% 0.05)
paper_dgp <- as.character(args$paper_dgp %||% "fMA1")
paper_T <- as.numeric(args$paper_T %||% 2)
paper_break_frac <- as.numeric(args$paper_break_frac %||% 0.8)
paper_m <- as.numeric(args$paper_m %||% 1000)
output_tier <- sanitize_tag(args$output_tier %||% "final_results")
output_root <- file.path(ROOT, "outputs", output_tier)

summary_file <- file.path(output_root, "summary", paste0("contam_training_summary_", scenario_target, ".csv"))
if (!file.exists(summary_file)) stop("Missing summary file: ", summary_file)

paper_root <- normalizePath(file.path(ROOT, ".."), winslash = "/", mustWork = TRUE)
fig_dir <- file.path(paper_root, "figs")
asset_dir <- file.path(output_root, "paper_assets")
ensure_dir(fig_dir)
ensure_dir(asset_dir)

sumdf <- read.csv(summary_file, stringsAsFactors = FALSE)
if (nrow(sumdf) == 0L) stop("Summary file is empty: ", summary_file)

for (nm in c("m", "T", "delta", "break_frac", "s_star", "alpha", "train_b",
             "train_break_frac", "train_drift_start_frac", "gamma", "bandwidth_h",
             "reject_rate", "post_break_detect_rate", "pre_break_stop_rate",
             "conditional_delay", "n_rep")) {
  if (nm %in% names(sumdf)) sumdf[[nm]] <- suppressWarnings(as.numeric(sumdf[[nm]]))
}
for (nm in c("dgp_type", "scenario", "train_contam", "contam_shape", "family",
             "standardizer", "detector", "type", "weight_name", "omega_name",
             "hset_name", "scale_weight_name", "method_group", "method_label",
             "method_id")) {
  if (nm %in% names(sumdf)) sumdf[[nm]] <- as.character(sumdf[[nm]])
}

target_df <- subset(
  sumdf,
  scenario == scenario_target &
    abs(alpha - alpha_target) < 1e-12 &
    abs(m - paper_m) < 1e-12 &
    abs(T - paper_T) < 1e-12 &
    abs(break_frac - paper_break_frac) < 1e-12
)
if (nrow(target_df) == 0L) {
  stop("No targeted rows found for scenario=", scenario_target,
       ", alpha=", alpha_target,
       ", m=", paper_m,
       ", T=", paper_T,
       ", break_frac=", paper_break_frac)
}

dgp_display <- function(x) {
  if (identical(x, "fMA1")) return("fMA(1)")
  x
}

main_ks_order <- c(
  "RSMS KS (gamma=0)",
  "RSMS KS (gamma=0.15)",
  "SSMS KS (gamma=0)",
  "SSMS KS (gamma=0.15)",
  "HAC KS (gamma=0)",
  "HAC KS (gamma=0.15)"
)

cvm_weight_blocks <- list(
  Uniform = c("RSMS weighted CvM [U]", "SSMS weighted CvM [U]", "HAC weighted CvM [U]"),
  Early   = c("RSMS weighted CvM [Early]", "SSMS weighted CvM [Early]", "HAC weighted CvM [Early]"),
  Mid     = c("RSMS weighted CvM [Mid]", "SSMS weighted CvM [Mid]", "HAC weighted CvM [Mid]"),
  Late    = c("RSMS weighted CvM [Late]", "SSMS weighted CvM [Late]", "HAC weighted CvM [Late]")
)

alt_cusum_order <- c(
  "RSMS WeightedCUSUM (gamma=0, InvSqrt)",
  "RSMS WeightedCUSUM (gamma=0.15, InvSqrt)",
  "SSMS WeightedCUSUM (gamma=0, InvSqrt)",
  "SSMS WeightedCUSUM (gamma=0.15, InvSqrt)",
  "HAC WeightedCUSUM (gamma=0, InvSqrt)",
  "HAC WeightedCUSUM (gamma=0.15, InvSqrt)",
  "RSMS PageCUSUM (gamma=0)",
  "RSMS PageCUSUM (gamma=0.15)",
  "SSMS PageCUSUM (gamma=0)",
  "SSMS PageCUSUM (gamma=0.15)",
  "HAC PageCUSUM (gamma=0)",
  "HAC PageCUSUM (gamma=0.15)"
)

alt_multiscale_order <- c(
  "RSMS MultiscaleMOSUM (H050_100_200, Equal)",
  "SSMS MultiscaleMOSUM (H050_100_200, Equal)",
  "HAC MultiscaleMOSUM (H050_100_200, Equal)"
)

style_for_method <- function(label) {
  if (grepl("^RSMS", label) && grepl("gamma=0.15", label)) return(list(col = "black", lty = 5, lwd = 2.2, pch = 1))
  if (grepl("^RSMS", label)) return(list(col = "black", lty = 1, lwd = 2.2, pch = 16))
  if (grepl("^SSMS", label) && grepl("gamma=0.15", label)) return(list(col = "gray35", lty = 4, lwd = 2.0, pch = 0))
  if (grepl("^SSMS", label)) return(list(col = "gray35", lty = 2, lwd = 2.0, pch = 15))
  if (grepl("^HAC", label) && grepl("gamma=0.15", label)) return(list(col = "gray60", lty = 6, lwd = 1.9, pch = 2))
  if (grepl("^HAC", label)) return(list(col = "gray60", lty = 3, lwd = 1.9, pch = 17))
  list(col = "black", lty = 1, lwd = 2, pch = 16)
}

style_for_multiscale <- function(label) {
  if (grepl("^RSMS", label)) return(list(col = "black", lty = 1, lwd = 2.3, pch = 16))
  if (grepl("^SSMS", label)) return(list(col = "gray35", lty = 2, lwd = 2.0, pch = 15))
  if (grepl("^HAC", label)) return(list(col = "gray60", lty = 3, lwd = 1.9, pch = 17))
  list(col = "black", lty = 1, lwd = 2, pch = 16)
}

legend_label_for_method <- function(label) {
  std <- if (grepl("^RSMS", label)) "RSMS" else if (grepl("^SSMS", label)) "SSMS" else if (grepl("^HAC", label)) "HAC" else label
  if (grepl("gamma=0.15", label)) return(sprintf("%s, gamma=0.15", std))
  if (grepl("gamma=0", label)) return(sprintf("%s, gamma=0", std))
  std
}

fmt_prob <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.1f", 100 * x))
}

fmt_delay <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.1f", x))
}

order_b_train <- function(x) sort(unique(x[is.finite(x)]))
order_contam <- function(x) intersect(c("late_break", "drift"), unique(x))

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}

lookup_metric <- function(df, method_name, contam_type, train_b_value, metric, dgp_value = NULL) {
  dd <- subset(df, method_label == method_name & train_contam == contam_type & abs(train_b - train_b_value) < 1e-12)
  if (!is.null(dgp_value)) dd <- subset(dd, dgp_type == dgp_value)
  if (nrow(dd) == 0L) return(NA_real_)
  dd[[metric]][1L]
}

plot_metric_family <- function(df, dgp_value, method_order, metric, outfile, title_text, subtitle_text = NULL, style_fun = style_for_method) {
  contam_levels <- order_contam(df$train_contam)
  b_vals <- order_b_train(df$train_b)
  png(outfile, width = 1600, height = 900, res = 150)
  par(mfrow = c(1, length(contam_levels)), mar = c(4.2, 4.4, 3.1, 1.1), oma = c(0, 0, 5.2, 0), xpd = NA)
  ylim <- c(0, 1)
  for (ct in contam_levels) {
    sub_df <- subset(df, dgp_type == dgp_value & train_contam == ct & method_label %in% method_order)
    plot(range(b_vals), ylim, type = "n", xlab = expression(b[train]), ylab = "", axes = FALSE,
         main = if (ct == "late_break") "Late-break contamination" else "Drift contamination")
    axis(1)
    axis(2, las = 1)
    box()
    grid(nx = NA, ny = NULL, col = "gray88", lty = 1)
    mtext(title_text, side = 2, line = 2.6, cex = 0.95)
    for (lab in method_order) {
      dd <- subset(sub_df, method_label == lab)
      if (nrow(dd) == 0L) next
      dd <- dd[order(dd$train_b), , drop = FALSE]
      sty <- style_fun(lab)
      lines(dd$train_b, dd[[metric]], col = sty$col, lty = sty$lty, lwd = sty$lwd)
      points(dd$train_b, dd[[metric]], col = sty$col, pch = sty$pch, cex = 0.9, bg = "white")
    }
  }
  legend("bottom", inset = c(0, -0.18), bty = "n", cex = 0.72,
         legend = vapply(method_order, legend_label_for_method, character(1)),
         col = vapply(method_order, function(z) style_fun(z)$col, character(1)),
         lty = vapply(method_order, function(z) style_fun(z)$lty, numeric(1)),
         lwd = vapply(method_order, function(z) style_fun(z)$lwd, numeric(1)),
         pch = vapply(method_order, function(z) style_fun(z)$pch, numeric(1)),
         ncol = 3)
  outer_title <- sprintf("%s (%s)", subtitle_text %||% title_text, dgp_display(dgp_value))
  mtext(outer_title, outer = TRUE, cex = 1.05, font = 2)
  dev.off()
}

make_two_block_table <- function(df, method_order, metric, caption, label, dgp_value, note = NULL) {
  b_vals <- order_b_train(df$train_b)
  nb <- length(b_vals)
  block_header <- paste(sprintf("& \\multicolumn{%d}{c}{%s}", nb, c("Late-break contamination", "Drift contamination")), collapse = " ")
  col_spec <- paste0("l", paste(rep("r", length(b_vals) * 2L), collapse = ""))

  make_rows <- function(labels) {
    rows <- character(0)
    for (lab in labels) {
      vals_late <- fmt_prob(vapply(b_vals, function(bb) lookup_metric(df, lab, "late_break", bb, metric, dgp_value), numeric(1)))
      vals_drift <- fmt_prob(vapply(b_vals, function(bb) lookup_metric(df, lab, "drift", bb, metric, dgp_value), numeric(1)))
      rows <- c(rows, paste0(paste(c(latex_escape(lab), vals_late, vals_drift), collapse = " & "), " \\\\"))
    }
    rows
  }

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{0.92}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{adjustbox}{max width=\\textwidth}",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    paste("Method", block_header, "\\\\"),
    sprintf("\\cmidrule(lr){2-%d} \\cmidrule(lr){%d-%d}", 1 + length(b_vals), 2 + length(b_vals), 1 + 2 * length(b_vals)),
    paste0(paste(c("Method", sprintf("$b_{\\mathrm{train}}=%s$", format(b_vals, trim = TRUE)), sprintf("$b_{\\mathrm{train}}=%s$", format(b_vals, trim = TRUE))), collapse = " & "), " \\\\"),
    "\\midrule"
  )
  lines <- c(lines, make_rows(method_order))
  lines <- c(lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}"
  )
  if (!is.null(note) && nzchar(note)) {
    lines <- c(lines, "\\begin{flushleft}", "\\footnotesize", note, "\\end{flushleft}")
  }
  lines <- c(lines, "\\end{table}")
  paste(lines, collapse = "\n")
}

make_weighted_cvm_table <- function(df, metric, caption, label, dgp_value, note = NULL) {
  b_vals <- order_b_train(df$train_b)
  nb <- length(b_vals)
  col_spec <- paste0("ll", paste(rep("r", length(b_vals) * 2L), collapse = ""))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.8pt}",
    "\\renewcommand{\\arraystretch}{0.92}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{adjustbox}{max width=\\textwidth}",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    sprintf("Weight & Method & \\multicolumn{%d}{c}{Late-break contamination} & \\multicolumn{%d}{c}{Drift contamination} \\\\", nb, nb),
    sprintf("\\cmidrule(lr){3-%d} \\cmidrule(lr){%d-%d}", 2 + nb, 3 + nb, 2 + 2 * nb),
    paste0(paste(c("Weight", "Method", sprintf("$b_{\\mathrm{train}}=%s$", format(b_vals, trim = TRUE)), sprintf("$b_{\\mathrm{train}}=%s$", format(b_vals, trim = TRUE))), collapse = " & "), " \\\\"),
    "\\midrule"
  )
  weight_labels <- names(cvm_weight_blocks)
  for (w in weight_labels) {
    method_block <- cvm_weight_blocks[[w]]
    for (j in seq_along(method_block)) {
      lab <- method_block[j]
      vals_late <- fmt_prob(vapply(b_vals, function(bb) lookup_metric(df, lab, "late_break", bb, metric, dgp_value), numeric(1)))
      vals_drift <- fmt_prob(vapply(b_vals, function(bb) lookup_metric(df, lab, "drift", bb, metric, dgp_value), numeric(1)))
      first_cell <- if (j == 1L) sprintf("\\multirow{3}{*}{%s}", w) else ""
      lines <- c(lines, paste0(paste(c(first_cell, sub(" weighted CvM \\[[^]]+\\]", "", latex_escape(lab)), vals_late, vals_drift), collapse = " & "), " \\\\"))
    }
    if (!identical(w, tail(weight_labels, 1L))) lines <- c(lines, "\\addlinespace[2pt]")
  }
  lines <- c(lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}"
  )
  if (!is.null(note) && nzchar(note)) {
    lines <- c(lines, "\\begin{flushleft}", "\\footnotesize", note, "\\end{flushleft}")
  }
  lines <- c(lines, "\\end{table}")
  paste(lines, collapse = "\n")
}

paper_df <- subset(target_df, dgp_type == paper_dgp)
if (nrow(paper_df) == 0L) stop("No paper-target rows for DGP=", paper_dgp)

main_ks_postbreak_fig <- file.path(fig_dir, "contam_main_ks_postbreak_fma1.png")
main_ks_prebreak_fig <- file.path(fig_dir, "contam_main_ks_prebreak_fma1.png")
alt_weighted_postbreak_fig <- file.path(fig_dir, "contam_alt_weighted_cusum_postbreak_fma1.png")
alt_weighted_prebreak_fig <- file.path(fig_dir, "contam_alt_weighted_cusum_prebreak_fma1.png")
alt_multiscale_postbreak_fig <- file.path(fig_dir, "contam_alt_multiscale_postbreak_fma1.png")
alt_multiscale_prebreak_fig <- file.path(fig_dir, "contam_alt_multiscale_prebreak_fma1.png")

plot_metric_family(paper_df, paper_dgp, main_ks_order, "post_break_detect_rate", main_ks_postbreak_fig,
                   title_text = "Post-break detection probability", subtitle_text = "Main KS robustness")
plot_metric_family(paper_df, paper_dgp, main_ks_order, "pre_break_stop_rate", main_ks_prebreak_fig,
                   title_text = "Pre-break stopping probability", subtitle_text = "Main KS robustness")
plot_metric_family(paper_df, paper_dgp, alt_cusum_order[1:6], "post_break_detect_rate", alt_weighted_postbreak_fig,
                   title_text = "Post-break detection probability", subtitle_text = "Weighted CUSUM robustness")
plot_metric_family(paper_df, paper_dgp, alt_cusum_order[1:6], "pre_break_stop_rate", alt_weighted_prebreak_fig,
                   title_text = "Pre-break stopping probability", subtitle_text = "Weighted CUSUM robustness")
plot_metric_family(paper_df, paper_dgp, alt_multiscale_order, "post_break_detect_rate", alt_multiscale_postbreak_fig,
                   title_text = "Post-break detection probability", subtitle_text = "Multiscale MOSUM robustness",
                   style_fun = style_for_multiscale)
plot_metric_family(paper_df, paper_dgp, alt_multiscale_order, "pre_break_stop_rate", alt_multiscale_prebreak_fig,
                   title_text = "Pre-break stopping probability", subtitle_text = "Multiscale MOSUM robustness",
                   style_fun = style_for_multiscale)

main_ks_table_tex <- make_two_block_table(
  paper_df,
  method_order = main_ks_order,
  metric = "post_break_detect_rate",
  caption = "Effect of mild training contamination on post-break detection for the main KS monitor family under the representative dependent design. Configuration: level-shift setting, fMA(1) errors, $m=1000$, $T=2$, $s^{\\star}/(mT)=0.8$, nominal level $\\alpha=0.05$, and Monte Carlo critical values from the clean-training calibration. Entries report $100\\times P(s^{\\star}<\\tau_m\\le mT)$.",
  label = "tab:contam-main-ks",
  dgp_value = paper_dgp,
  note = "Note: The left block injects a late break inside the training sample at $\\lfloor 0.8m\\rfloor$, while the right block imposes a linear training drift whose endpoint level persists into Phase II."
)

main_cvm_table_tex <- make_weighted_cvm_table(
  paper_df,
  metric = "post_break_detect_rate",
  caption = "Effect of mild training contamination on post-break detection for the weighted-CvM monitor family under the representative dependent design. Configuration: level-shift setting, fMA(1) errors, $m=1000$, $T=2$, $s^{\\star}/(mT)=0.8$, nominal level $\\alpha=0.05$, and clean-training critical values.",
  label = "tab:contam-main-cvm",
  dgp_value = paper_dgp,
  note = "Note: Within each weight block, the rows list RSMS, SSMS, and HAC. Entries report $100\\times P(s^{\\star}<\\tau_m\\le mT)$."
)

alt_cusum_table_tex <- make_two_block_table(
  paper_df,
  method_order = alt_cusum_order,
  metric = "post_break_detect_rate",
  caption = "Effect of mild training contamination on post-break detection for the retained CUSUM-style alternative detectors under the representative dependent design. Configuration: level-shift setting, fMA(1) errors, $m=1000$, $T=2$, $s^{\\star}/(mT)=0.8$, and $\\alpha=0.05$.",
  label = "tab:contam-alt-cusum",
  dgp_value = paper_dgp,
  note = "Note: The weighted-CUSUM rows use the fixed segment-length weight $\\omega(\\ell)=\\ell^{-1/2}$. The Page-CUSUM rows are retained as a restart benchmark but are generally weaker than weighted CUSUM in the paper's narrative."
)

alt_multiscale_table_tex <- make_two_block_table(
  paper_df,
  method_order = alt_multiscale_order,
  metric = "post_break_detect_rate",
  caption = "Effect of mild training contamination on post-break detection for the retained multiscale MOSUM benchmark under the representative dependent design. Configuration: level-shift setting, fMA(1) errors, $m=1000$, $T=2$, $s^{\\star}/(mT)=0.8$, and $\\alpha=0.05$.",
  label = "tab:contam-alt-multiscale",
  dgp_value = paper_dgp,
  note = "Note: The multiscale scan uses $\\mathcal H=\\{0.05,0.10,0.20\\}$ with equal scale weights."
)

atomic_write_lines(strsplit(main_ks_table_tex, "\n", fixed = TRUE)[[1L]], file.path(asset_dir, "contam_main_ks_table.tex"))
atomic_write_lines(strsplit(main_cvm_table_tex, "\n", fixed = TRUE)[[1L]], file.path(asset_dir, "contam_main_cvm_table.tex"))
atomic_write_lines(strsplit(alt_cusum_table_tex, "\n", fixed = TRUE)[[1L]], file.path(asset_dir, "contam_alt_cusum_table.tex"))
atomic_write_lines(strsplit(alt_multiscale_table_tex, "\n", fixed = TRUE)[[1L]], file.path(asset_dir, "contam_alt_multiscale_table.tex"))

avg_metric <- function(df, rows, metric, train_b_max = 0.005) {
  dd <- subset(df, method_label %in% rows & train_b <= train_b_max)
  if (nrow(dd) == 0L) return(NA_real_)
  mean(dd[[metric]], na.rm = TRUE)
}

rsms_ks_mild <- avg_metric(paper_df, c("RSMS KS (gamma=0)", "RSMS KS (gamma=0.15)"), "post_break_detect_rate", 0.005)
ssms_ks_mild <- avg_metric(paper_df, c("SSMS KS (gamma=0)", "SSMS KS (gamma=0.15)"), "post_break_detect_rate", 0.005)
hac_ks_mild <- avg_metric(paper_df, c("HAC KS (gamma=0)", "HAC KS (gamma=0.15)"), "post_break_detect_rate", 0.005)
rsms_ks_all_dgp_mild <- avg_metric(target_df, c("RSMS KS (gamma=0)", "RSMS KS (gamma=0.15)"), "post_break_detect_rate", 0.005)
ssms_ks_all_dgp_mild <- avg_metric(target_df, c("SSMS KS (gamma=0)", "SSMS KS (gamma=0.15)"), "post_break_detect_rate", 0.005)
hac_ks_all_dgp_mild <- avg_metric(target_df, c("HAC KS (gamma=0)", "HAC KS (gamma=0.15)"), "post_break_detect_rate", 0.005)

rsms_wc_mild <- avg_metric(paper_df, c("RSMS WeightedCUSUM (gamma=0, InvSqrt)", "RSMS WeightedCUSUM (gamma=0.15, InvSqrt)"), "post_break_detect_rate", 0.005)
rsms_page_mild <- avg_metric(paper_df, c("RSMS PageCUSUM (gamma=0)", "RSMS PageCUSUM (gamma=0.15)"), "post_break_detect_rate", 0.005)
rsms_ms_mild <- avg_metric(paper_df, c("RSMS MultiscaleMOSUM (H050_100_200, Equal)"), "post_break_detect_rate", 0.005)

notes <- c(
  sprintf("Targeted robustness design: scenario=%s, dgp=%s, m=%s, T=%s, break_frac=%.2f, alpha=%.2f.",
          scenario_target, dgp_display(paper_dgp), paper_m, paper_T, paper_break_frac, alpha_target),
  sprintf("Main KS mild-contamination averages (train_b<=0.005): RSMS=%.1f%%, SSMS=%.1f%%, HAC=%.1f%% post-break detection.",
          100 * rsms_ks_mild, 100 * ssms_ks_mild, 100 * hac_ks_mild),
  sprintf("Cross-DGP KS mild-contamination averages (train_b<=0.005): RSMS=%.1f%%, SSMS=%.1f%%, HAC=%.1f%% post-break detection.",
          100 * rsms_ks_all_dgp_mild, 100 * ssms_ks_all_dgp_mild, 100 * hac_ks_all_dgp_mild),
  sprintf("Alternative mild-contamination averages (train_b<=0.005): RSMS weighted CUSUM=%.1f%%, RSMS Page-CUSUM=%.1f%%, RSMS multiscale MOSUM=%.1f%% post-break detection.",
          100 * rsms_wc_mild, 100 * rsms_page_mild, 100 * rsms_ms_mild),
  sprintf("Selected figures and tables are written here; CSV-level summaries remain in outputs/%s/summary.", output_tier)
)
atomic_write_lines(notes, file.path(asset_dir, "contam_paper_notes.txt"))

message("Wrote contamination-study tables and figures to: ", normalize_path2(asset_dir, mustWork = FALSE))
message("Updated figures in: ", normalize_path2(fig_dir, mustWork = FALSE))
