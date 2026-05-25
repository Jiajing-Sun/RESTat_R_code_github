#!/usr/bin/env Rscript

# ==============================================================
# build_hac_tuning_sensitivity.R
# Supplementary study of HAC kernel and bandwidth sensitivity
# for the main KS statistic under the fMA(1) setting.
# ==============================================================

decode_rscript_path <- function(x) gsub("~\\+~", " ", x)
file_args <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- normalizePath(decode_rscript_path(sub("^--file=", "", file_args[1L])), mustWork = FALSE)
script_dir <- dirname(script_path)
main_stats_root <- dirname(script_dir)
bundle_root <- dirname(main_stats_root)
project_root <- dirname(bundle_root)
workspace_root <- dirname(project_root)
cv_root <- main_stats_root
rmod_dir <- file.path(main_stats_root, "R")

source(file.path(rmod_dir, "dependencies.R"), local = TRUE)
source(file.path(rmod_dir, "project_paths.R"), local = TRUE)
source(file.path(rmod_dir, "utils.R"), local = TRUE)
source(file.path(rmod_dir, "simulation_settings.R"), local = TRUE)
source(file.path(rmod_dir, "genData.R"), local = TRUE)
source(file.path(rmod_dir, "scenarios.R"), local = TRUE)
source(file.path(rmod_dir, "fpca_pipeline.R"), local = TRUE)
source(file.path(rmod_dir, "critical_values_lookup.R"), local = TRUE)
source(file.path(rmod_dir, "detectors_main.R"), local = TRUE)

ensure_packages(c("fdapace", "fda", "data.table", "ggplot2"), install_if_missing = FALSE)

args <- parse_named_args()
nrep <- as.integer(args$nrep %||% 1000L)
ncores <- as.integer(args$ncores %||% max(1L, parallel::detectCores(logical = TRUE) - 1L))
seed0 <- as.integer(args$seed %||% 20260425L)

bundle_output_dir <- file.path(main_stats_root, "outputs")
bundle_figure_dir <- file.path(main_stats_root, "figures")
submission_figure_dir <- bundle_figure_dir
submission_generated_dir <- file.path(main_stats_root, "generated")

ensure_dir(bundle_output_dir)
ensure_dir(bundle_figure_dir)
ensure_dir(submission_figure_dir)
ensure_dir(submission_generated_dir)

summary_csv <- file.path(bundle_output_dir, "hac_tuning_sensitivity_summary.csv")
raw_csv <- file.path(bundle_output_dir, "hac_tuning_sensitivity_raw.csv")
figure_bundle <- file.path(bundle_figure_dir, "hac_tuning_sensitivity_tradeoff.png")
figure_submission <- file.path(submission_figure_dir, "hac_tuning_sensitivity_tradeoff.png")
table_tex <- file.path(submission_generated_dir, "hac_tuning_sensitivity_table.tex")

hac_kernel_weight <- function(ell, bandwidth, kernel) {
  kernel <- match.arg(kernel, c("Bartlett", "Parzen", "QS"))
  x <- ell / bandwidth
  if (kernel == "Bartlett") {
    if (x > 1) return(0)
    return(1 - x)
  }
  if (kernel == "Parzen") {
    if (x > 1) return(0)
    if (x <= 0.5) return(1 - 6 * x^2 + 6 * x^3)
    return(2 * (1 - x)^3)
  }
  if (x == 0) return(1)
  a <- 6 * pi * x / 5
  25 / (12 * pi^2 * x^2) * (sin(a) / a - cos(a))
}

hac_long_run_covariance <- function(Z, bandwidth = NULL, kernel = "Bartlett", center = TRUE, ridge = 1e-10) {
  Z <- safe_matrix(Z)
  m <- nrow(Z)
  if (m < 2L) return(diag(rep(1, ncol(Z))))
  if (center) Z <- sweep(Z, 2, colMeans(Z), FUN = "-")
  L <- if (is.null(bandwidth)) hac_bandwidth_default(m) else as.integer(max(1L, round(bandwidth[1L])))
  L <- min(L, m - 1L)
  S <- crossprod(Z) / m
  lag_grid <- if (identical(kernel, "QS")) seq_len(m - 1L) else seq_len(L)
  for (ell in lag_grid) {
    w <- hac_kernel_weight(ell, L, kernel = kernel)
    if (!is.finite(w) || abs(w) < 1e-12) next
    Gell <- crossprod(Z[(ell + 1L):m, , drop = FALSE], Z[1L:(m - ell), , drop = FALSE]) / m
    S <- S + w * (Gell + t(Gell))
  }
  (S + t(S)) / 2 + ridge * diag(ncol(Z))
}

evaluate_path_from_context <- function(ctx, q_used, gamma, critical_value) {
  s <- ctx$process$s
  k_index <- seq_len(ctx$process$k_max)
  quad <- quad_metric_rows(ctx$process$U_full[, 1:q_used, drop = FALSE], ctx)
  path <- quad / g_gamma_sq(s, gamma)
  data.frame(
    statistic = max(path),
    reject = max(path) > critical_value,
    first_rejection = first_crossing(path, critical_value, k_index = k_index),
    stringsAsFactors = FALSE
  )
}

prepare_hac_context_custom <- function(scores_train, scores_monitor, kernel, bandwidth, ridge = 1e-10) {
  centered <- center_train_monitor(scores_train, scores_monitor)
  proc <- build_score_processes(centered$train, centered$monitor)
  Sigma <- hac_long_run_covariance(centered$train, bandwidth = bandwidth, kernel = kernel, center = FALSE, ridge = ridge)
  list(
    standardizer = "HAC",
    process = proc,
    metric = list(kind = "matrix", inv = safe_solve(Sigma, ridge = ridge))
  )
}

simulate_one_replicate <- function(rep_id, scenario_name, delta, s_star, gamma_vec, hac_specs, cv_main, seed) {
  set.seed(seed)
  sim_m <- 500L
  sim_T <- 2
  alpha <- 0.05
  ridge <- 1e-10

  total_n <- sim_m + as.integer(round(sim_m * sim_T))
  fd_data <- generate_base_fd("fMA1", total_n = total_n)
  fd_data <- apply_change_to_fd(
    fd_data,
    scenario = scenario_name,
    m = sim_m,
    T = sim_T,
    delta = delta,
    s_star = s_star,
    basis_k = localized_change_basis_index_default()
  )

  fpca <- fpca_project_scores(
    fd_data,
    m = sim_m,
    q_cap = 30L,
    fve_threshold = 0.95,
    fixed_q = NA_integer_,
    n_grid = 301L
  )
  q_used <- fpca$q_used
  whitened <- whiten_scores(fpca$scores_train, fpca$scores_monitor, ridge = 1e-8)

  ctx_ss <- prepare_standardizer_context(fpca$scores_train, fpca$scores_monitor, "SSMS", ridge = ridge, range_floor = 1e-8)
  ctx_rs <- prepare_standardizer_context(whitened$scores_train, whitened$scores_monitor, "RSMS", ridge = ridge, range_floor = 1e-8)

  out <- list()
  idx <- 1L
  for (gamma in gamma_vec) {
    cv_ss <- lookup_main_critical_value(cv_main, standardizer = "SSMS", type = "KS", T = sim_T, q = q_used, alpha = alpha, gamma = gamma)
    cv_rs <- lookup_main_critical_value(cv_main, standardizer = "RSMS", type = "KS", T = sim_T, q = q_used, alpha = alpha, gamma = gamma)
    cv_hac <- lookup_main_critical_value(cv_main, standardizer = "HAC", type = "KS", T = sim_T, q = q_used, alpha = alpha, gamma = gamma)

    rs_eval <- evaluate_path_from_context(ctx_rs, q_used = q_used, gamma = gamma, critical_value = cv_rs)
    out[[idx]] <- data.frame(
      rep = rep_id,
      scenario = scenario_name,
      delta = delta,
      s_star = s_star,
      gamma = gamma,
      method_id = sprintf("rsms_ks_gamma_%s", gsub("\\.", "_", format(gamma, trim = TRUE))),
      paper_label = "RSMS KS",
      standardizer = "RSMS",
      kernel = NA_character_,
      bandwidth = NA_integer_,
      bandwidth_label = "Anchor",
      statistic = rs_eval$statistic,
      reject = rs_eval$reject,
      first_rejection = rs_eval$first_rejection,
      q_used = q_used,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L

    ss_eval <- evaluate_path_from_context(ctx_ss, q_used = q_used, gamma = gamma, critical_value = cv_ss)
    out[[idx]] <- data.frame(
      rep = rep_id,
      scenario = scenario_name,
      delta = delta,
      s_star = s_star,
      gamma = gamma,
      method_id = sprintf("ssms_ks_gamma_%s", gsub("\\.", "_", format(gamma, trim = TRUE))),
      paper_label = "SSMS KS",
      standardizer = "SSMS",
      kernel = NA_character_,
      bandwidth = NA_integer_,
      bandwidth_label = "Anchor",
      statistic = ss_eval$statistic,
      reject = ss_eval$reject,
      first_rejection = ss_eval$first_rejection,
      q_used = q_used,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L

    for (j in seq_len(nrow(hac_specs))) {
      spec <- hac_specs[j, , drop = FALSE]
      ctx_h <- prepare_hac_context_custom(
        fpca$scores_train,
        fpca$scores_monitor,
        kernel = spec$kernel,
        bandwidth = spec$bandwidth,
        ridge = ridge
      )
      h_eval <- evaluate_path_from_context(ctx_h, q_used = q_used, gamma = gamma, critical_value = cv_hac)
      out[[idx]] <- data.frame(
        rep = rep_id,
        scenario = scenario_name,
        delta = delta,
        s_star = s_star,
        gamma = gamma,
        method_id = sprintf("hac_%s_L%d_gamma_%s", tolower(spec$kernel), spec$bandwidth, gsub("\\.", "_", format(gamma, trim = TRUE))),
        paper_label = sprintf("HAC %s (L=%d%s)", spec$kernel, spec$bandwidth, ifelse(spec$is_default, ", default", "")),
        standardizer = "HAC",
        kernel = spec$kernel,
        bandwidth = spec$bandwidth,
        bandwidth_label = sprintf("L=%d", spec$bandwidth),
        statistic = h_eval$statistic,
        reject = h_eval$reject,
        first_rejection = h_eval$first_rejection,
        q_used = q_used,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

fmt_pct <- function(x) sprintf("%.1f", 100 * x)
fmt_num <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))

write_table_tex <- function(summary_dt, path) {
  order_levels <- c("RSMS KS", "SSMS KS")
  ordered_methods <- c(
    order_levels,
    summary_dt[standardizer == "HAC" & gamma == 0, unique(paper_label)]
  )
  panel_titles <- c("0" = "\\textbf{Panel A:} $\\gamma=0$", "0.15" = "\\textbf{Panel B:} $\\gamma=0.15$")

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{0.95}",
    "\\caption{HAC kernel and bandwidth sensitivity for the main KS statistic under the fMA(1) setting.}",
    "\\label{tab:hac-tuning-sensitivity}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    "Method & Size & Level Raw & Level SAP & Level ADD & Smooth Raw & Smooth SAP & Smooth ADD \\\\",
    "\\midrule"
  )

  for (g in c(0, 0.15)) {
    sub <- summary_dt[gamma == g]
    sub$order_id <- match(sub$paper_label, ordered_methods)
    sub <- sub[order(sub$order_id), ]
    lines <- c(lines, sprintf("\\multicolumn{8}{l}{%s} \\\\", panel_titles[as.character(g)]), "\\midrule")
    for (ii in seq_len(nrow(sub))) {
      lines <- c(
        lines,
        sprintf(
          "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
          sub$paper_label[ii],
          fmt_pct(sub$size[ii]),
          fmt_pct(sub$level_raw[ii]),
          fmt_pct(sub$level_sap[ii]),
          fmt_num(sub$level_add[ii]),
          fmt_pct(sub$smooth_raw[ii]),
          fmt_pct(sub$smooth_sap[ii]),
          fmt_num(sub$smooth_add[ii])
        )
      )
    }
    if (!identical(g, 0.15)) lines <- c(lines, "\\midrule")
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\begin{flushleft}",
    sprintf(
      "Note: The study adopts the fMA(1) setting used in the paper: $m=500$, $T=2$, $s^{\\star}=200$, and 95\\%% FVE FPCA compression. Empirical size is evaluated at nominal 5\\%% under $\\Delta=0$. The level-shift and smooth-change columns use $\\Delta=0.017$. Raw power uses the asymptotic 5\\%% critical value. Size-adjusted power (SAP) uses each method's empirical 95th percentile under the null. Average detection delay (ADD) is the mean post-break stopping delay among asymptotic rejections. The HAC rows use kernel-specific long-run covariance estimators with bandwidth $L$; the default Bartlett bandwidth for $m=500$ is $L=%d$.",
      hac_bandwidth_default(500L)
    ),
    "\\end{flushleft}",
    "\\end{table}"
  )
  write_lines_atomic(lines, path)
}

default_bw <- hac_bandwidth_default(500L)
hac_specs <- data.table::CJ(
  kernel = c("Bartlett", "Parzen", "QS"),
  bandwidth = c(max(1L, floor(default_bw / 2)), default_bw, min(500L - 1L, 2L * default_bw))
)
hac_specs[, is_default := bandwidth == default_bw]

gamma_vec <- c(0, 0.15)
scenario_grid <- data.table::data.table(
  scenario = c("null", "level_shift", "smooth_change"),
  delta = c(0, 0.017, 0.017),
  s_star = c(200L, 200L, 200L)
)

cv_main <- load_main_critical_values(cv_root)
seed_stream <- seed0 + seq_len(nrep)
rep_grid <- expand.grid(rep = seq_len(nrep), scenario_idx = seq_len(nrow(scenario_grid)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

worker_fun <- function(ii) {
  rep_id <- rep_grid$rep[ii]
  ss <- scenario_grid[rep_grid$scenario_idx[ii]]
  simulate_one_replicate(
    rep_id = rep_id,
    scenario_name = ss$scenario[[1L]],
    delta = ss$delta[[1L]],
    s_star = ss$s_star[[1L]],
    gamma_vec = gamma_vec,
    hac_specs = hac_specs,
    cv_main = cv_main,
    seed = seed_stream[rep_id] + 10000L * rep_grid$scenario_idx[ii]
  )
}

message("Running focused HAC tuning study with nrep=", nrep, " and ncores=", ncores, " ...")

res_list <- if (.Platform$OS.type == "windows" || ncores <= 1L) {
  lapply(seq_len(nrow(rep_grid)), worker_fun)
} else {
  parallel::mclapply(seq_len(nrow(rep_grid)), worker_fun, mc.cores = ncores, mc.preschedule = TRUE)
}

raw_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
data.table::setDT(raw_dt)
write_csv_atomic(raw_dt, raw_csv, row.names = FALSE)

null_dt <- raw_dt[scenario == "null", .(
  size = mean(reject, na.rm = TRUE),
  empirical_cv = as.numeric(stats::quantile(statistic, probs = 0.95, type = 7, na.rm = TRUE))
), by = .(gamma, method_id, paper_label, standardizer, kernel, bandwidth, bandwidth_label)]

power_dt <- merge(
  raw_dt[scenario %in% c("level_shift", "smooth_change")],
  null_dt[, .(gamma, method_id, empirical_cv)],
  by = c("gamma", "method_id"),
  all.x = TRUE,
  sort = FALSE
)
power_dt[, size_adjusted_reject := statistic > empirical_cv]
power_dt[, delay := ifelse(reject & !is.na(first_rejection) & first_rejection >= s_star, first_rejection - s_star, NA_real_)]

power_sum <- power_dt[, .(
  raw_power = mean(reject, na.rm = TRUE),
  sap = mean(size_adjusted_reject, na.rm = TRUE),
  add = if (all(is.na(delay))) NA_real_ else mean(delay, na.rm = TRUE)
), by = .(gamma, method_id, paper_label, standardizer, kernel, bandwidth, bandwidth_label, scenario)]

summary_dt <- merge(
  null_dt,
  power_sum[scenario == "level_shift", .(
    gamma, method_id, level_raw = raw_power, level_sap = sap, level_add = add
  )],
  by = c("gamma", "method_id"),
  all.x = TRUE,
  sort = FALSE
)
summary_dt <- merge(
  summary_dt,
  power_sum[scenario == "smooth_change", .(
    gamma, method_id, smooth_raw = raw_power, smooth_sap = sap, smooth_add = add
  )],
  by = c("gamma", "method_id"),
  all.x = TRUE,
  sort = FALSE
)
summary_dt <- summary_dt[order(gamma, standardizer, kernel, bandwidth)]
write_csv_atomic(summary_dt, summary_csv, row.names = FALSE)
write_table_tex(summary_dt, table_tex)

plot_dt <- merge(
  power_sum[scenario == "level_shift", .(
    gamma, method_id, paper_label, standardizer, kernel, bandwidth, bandwidth_label,
    metric = "Raw Power", value = raw_power
  )],
  null_dt[, .(gamma, method_id, size)],
  by = c("gamma", "method_id"),
  all.x = TRUE,
  sort = FALSE
)
plot_dt <- data.table::rbindlist(list(
  plot_dt,
  merge(
    power_sum[scenario == "level_shift", .(
      gamma, method_id, paper_label, standardizer, kernel, bandwidth, bandwidth_label,
      metric = "SAP", value = sap
    )],
    null_dt[, .(gamma, method_id, size)],
    by = c("gamma", "method_id"),
    all.x = TRUE,
    sort = FALSE
  )
), use.names = TRUE, fill = TRUE)
plot_dt <- data.table::rbindlist(list(
  plot_dt,
  null_dt[, .(
    gamma, method_id, paper_label, standardizer, kernel, bandwidth, bandwidth_label,
    metric = "Size", value = size, size = size
  )]
), use.names = TRUE, fill = TRUE)

anchors <- plot_dt[standardizer %in% c("RSMS", "SSMS"), .(
  gamma, metric, paper_label, value
)]
hac_plot <- plot_dt[standardizer == "HAC"]
kernel_cols <- c(Bartlett = "#0072B2", Parzen = "#D55E00", QS = "#009E73")
anchor_cols <- c("RSMS KS" = "black", "SSMS KS" = "gray45")

p <- ggplot2::ggplot(hac_plot, ggplot2::aes(x = bandwidth, y = value, color = kernel, shape = kernel, group = kernel)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::geom_hline(
    data = anchors,
    ggplot2::aes(yintercept = value, linetype = paper_label),
    color = "gray20",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_grid(metric ~ gamma, labeller = ggplot2::labeller(
    gamma = c(`0` = expression(gamma == 0), `0.15` = expression(gamma == 0.15))
  ), scales = "free_y") +
  ggplot2::scale_color_manual(values = kernel_cols) +
  ggplot2::scale_shape_manual(values = c(Bartlett = 16, Parzen = 17, QS = 15)) +
  ggplot2::scale_linetype_manual(values = c("RSMS KS" = "dashed", "SSMS KS" = "dotted")) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(hac_specs$bandwidth))
  ) +
  ggplot2::labs(
    x = "HAC bandwidth L",
    y = NULL,
    title = "HAC Kernel-and-Bandwidth Sensitivity Analysis for the Main KS Monitor",
    subtitle = "Representative dependent design: fMA(1), m=500, T=2, s*=200, and Delta=0.017 for the level shift"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(size = 10.5),
    legend.text = ggplot2::element_text(size = 10),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    strip.background = ggplot2::element_rect(fill = "#DCEEFF", color = "gray60"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(title = "HAC kernel", nrow = 1, byrow = TRUE, order = 1),
    shape = ggplot2::guide_legend(title = "HAC kernel", nrow = 1, byrow = TRUE, order = 1),
    linetype = ggplot2::guide_legend(title = "Benchmark line", nrow = 1, byrow = TRUE, order = 2)
  )

ggplot2::ggsave(figure_bundle, p, width = 11.0, height = 9.0, dpi = 220)
file.copy(figure_bundle, figure_submission, overwrite = TRUE)

message("Saved focused HAC tuning summary to: ", summary_csv)
message("Saved focused HAC tuning table to: ", table_tex)
message("Saved focused HAC tuning figure to: ", figure_submission)
