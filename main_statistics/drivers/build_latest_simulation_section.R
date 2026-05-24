suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

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

root_dir_override <- Sys.getenv("STREAMINGCURVE_PAPER_ROOT", unset = "")
root_dir <- if (nzchar(root_dir_override)) {
  normalizePath(root_dir_override, winslash = "/", mustWork = FALSE)
} else {
  dirname(bootstrap_script_path())
}
sim_root <- file.path(root_dir, "fresh_streaming_curve_sim_codes_副本")
null_dir <- file.path(sim_root, "outputs", "null_raw")
power_root <- file.path(sim_root, "outputs", "power_raw")
paper_summary_dir <- file.path(sim_root, "outputs", "paper_summary")
fig_dir <- file.path(root_dir, "figs")
tex_path <- file.path(root_dir, "Online_Monitoring_via_Streaming_Curves.tex")

dir.create(paper_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

main_scenarios <- c("level_shift", "smooth_change")
appendix_main_scenarios <- c("abrupt_local_change")
all_scenarios <- c(main_scenarios, appendix_main_scenarios)
weight_levels <- c("U", "Early", "Mid", "Late")
standardizer_levels <- c("RSMS", "SSMS", "HAC")
dgp_levels <- c("BB", "fIID", "fMA1")
t_levels <- c("1", "2", "5", "10")
alt_detector_keep <- c("WeightedCUSUM", "PageCUSUM", "MultiscaleMOSUM")

marker_main_start <- "% AUTO-GENERATED MAIN SIMULATION SECTION START"
marker_main_end <- "% AUTO-GENERATED MAIN SIMULATION SECTION END"
marker_app_start <- "% AUTO-GENERATED APPENDIX SIMULATION SECTION START"
marker_app_end <- "% AUTO-GENERATED APPENDIX SIMULATION SECTION END"

selected_cols <- c(
  "family", "standardizer", "detector", "type", "gamma", "weight_name",
  "bandwidth_h", "hset_name", "scale_weight_name", "method_group",
  "method_label", "method_id", "alpha", "statistic", "reject",
  "first_rejection", "dgp_type", "scenario", "m", "T", "delta", "s_star",
  "error_flag"
)

scenario_pretty <- c(
  level_shift = "Level shift",
  smooth_change = "Smooth change",
  abrupt_local_change = "Localized change"
)
scenario_pretty_adj <- c(
  level_shift = "level-shift",
  smooth_change = "smooth-change",
  abrupt_local_change = "localized-change"
)
dgp_pretty <- c(BB = "BB", fIID = "fIID", fMA1 = "fMA(1)")
weight_pretty <- c(
  U = "Uniform weight $w_{\\mathrm U}$",
  Early = "Early-emphasis weight $w_{\\mathrm E}$",
  Mid = "Mid-emphasis weight $w_{\\mathrm M}$",
  Late = "Late-emphasis weight $w_{\\mathrm L}$"
)
weight_caption_pretty <- c(
  U = "the uniform weight $w_{\\mathrm U}$",
  Early = "the early-emphasis weight $w_{\\mathrm E}$",
  Mid = "the mid-emphasis weight $w_{\\mathrm M}$",
  Late = "the late-emphasis weight $w_{\\mathrm L}$"
)
weight_short_pretty <- c(
  U = "uniform",
  Early = "early-emphasis",
  Mid = "mid-emphasis",
  Late = "late-emphasis"
)
weight_table_pretty <- c(
  U = "Uniform",
  Early = "Early",
  Mid = "Mid",
  Late = "Late"
)
detector_pretty <- c(
  PageCUSUM = "Page-CUSUM",
  WeightedCUSUM = "weighted CUSUM",
  MultiscaleMOSUM = "multiscale MOSUM"
)

clean_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- trimws(toupper(as.character(x)))
  y %in% c("TRUE", "T", "1", "YES", "Y")
}

gamma_label <- function(x) {
  if (is.na(x)) return("")
  if (abs(x) < 1e-12) return("0")
  formatC(x, format = "f", digits = 2, drop0trailing = TRUE)
}

fmt_pct <- function(x) ifelse(is.finite(x), sprintf("%.1f", 100 * x), "--")
fmt_pp <- function(x) ifelse(is.finite(x), sprintf("%.1f", 100 * x), "--")
fmt_num <- function(x) ifelse(is.finite(x), sprintf("%.1f", x), "--")
fmt_delta <- function(x) formatC(x, format = "f", digits = 3, drop0trailing = TRUE)
t_numeric_sort <- function(x) as.character(sort(unique(as.numeric(as.character(x)))))

fig_rel <- function(name) file.path("figs", name)

main_method_label <- function(standardizer, type, gamma, weight_name) {
  if (identical(type, "KS")) {
    sprintf("%s KS, $\\gamma=%s$", standardizer, gamma_label(gamma))
  } else {
    sprintf("%s CvM-%s", standardizer, weight_name)
  }
}

alt_tuning_label <- function(detector, gamma, bandwidth_h, hset_name) {
  if (identical(detector, "PageCUSUM")) {
    return(sprintf("Page-CUSUM, $\\gamma=%s$", gamma_label(gamma)))
  }
  if (identical(detector, "WeightedCUSUM")) {
    return(sprintf("weighted CUSUM, $\\gamma=%s$", gamma_label(gamma)))
  }
  if (identical(detector, "MultiscaleMOSUM")) {
    if (!is.na(hset_name) && nzchar(hset_name)) {
      return(sprintf("Multiscale MOSUM, %s", gsub("_", "\\_", hset_name, fixed = TRUE)))
    }
    return("Multiscale MOSUM")
  }
  detector
}

alt_tuning_label_plot <- function(detector, gamma, bandwidth_h, hset_name) {
  if (identical(detector, "PageCUSUM")) {
    return(sprintf("Page-CUSUM, gamma=%s", gamma_label(gamma)))
  }
  if (identical(detector, "WeightedCUSUM")) {
    return(sprintf("weighted CUSUM, gamma=%s", gamma_label(gamma)))
  }
  if (identical(detector, "MultiscaleMOSUM")) {
    if (!is.na(hset_name) && nzchar(hset_name)) {
      return(sprintf("Multiscale MOSUM, %s", hset_name))
    }
    return("Multiscale MOSUM")
  }
  detector
}

build_method_meta <- function(summary_dt, family_tag) {
  meta <- unique(summary_dt[, .(
    family, standardizer, detector, type, gamma, weight_name,
    bandwidth_h, hset_name, scale_weight_name, method_id
  )])

  if (identical(family_tag, "Main")) {
    meta[, weight_name := fifelse(is.na(weight_name), "", weight_name)]
    meta[, paper_label := mapply(
      main_method_label, standardizer, type, gamma, weight_name,
      USE.NAMES = FALSE
    )]
    meta[, panel := fifelse(type == "KS", "KS monitors", "Weighted CvM monitors")]
    meta[, std_rank := match(standardizer, standardizer_levels)]
    meta[, weight_rank := fifelse(type == "CvM", match(weight_name, weight_levels), 0L)]
    meta[, gamma_rank := fifelse(is.na(gamma), 0L, fifelse(abs(gamma - 0.15) < 1e-12, 2L, 1L))]
    meta[, method_rank := fifelse(
      type == "KS",
      std_rank * 10L + gamma_rank,
      100L + std_rank * 10L + weight_rank
    )]
  } else {
    meta[, tuning_label := mapply(
      alt_tuning_label, detector, gamma, bandwidth_h, hset_name,
      USE.NAMES = FALSE
    )]
    meta[, tuning_label_plot := mapply(
      alt_tuning_label_plot, detector, gamma, bandwidth_h, hset_name,
      USE.NAMES = FALSE
    )]
    meta[, paper_label := sprintf("%s %s", standardizer, tuning_label)]
    meta[, panel := "Alternative detectors"]
    meta[, std_rank := match(standardizer, standardizer_levels)]
    meta[, detector_rank := fifelse(
      detector == "WeightedCUSUM",
      10L + fifelse(abs(gamma - 0.15) < 1e-12, 2L, 1L),
      fifelse(
        detector == "PageCUSUM",
        20L + fifelse(abs(gamma - 0.15) < 1e-12, 2L, 1L),
        30L
      )
    )]
    meta[, method_rank := detector_rank * 10L + std_rank]
  }

  setorder(meta, method_rank)
  meta
}

filter_main_dt <- function(dt) {
  dt[
    family == "Main" &
      abs(alpha - 0.05) < 1e-12 &
      !clean_bool(error_flag)
  ]
}

filter_alt_dt <- function(dt) {
  dt[
    family == "Benchmark" &
      detector %in% alt_detector_keep &
      abs(alpha - 0.05) < 1e-12 &
      !clean_bool(error_flag)
  ]
}

summarize_null_files <- function(paths, filter_fun) {
  rbindlist(lapply(paths, function(path) {
    dt <- fread(path, select = selected_cols, showProgress = FALSE)
    dt <- filter_fun(dt)
    if (nrow(dt) == 0L) return(NULL)
    dt[, reject := clean_bool(reject)]
    dt[, T := as.character(T)]
    dt[
      ,
      .(
        size = mean(reject),
        empirical_cv = as.numeric(quantile(statistic, probs = 0.95, type = 7, na.rm = TRUE)),
        n_rep = .N
      ),
      by = .(
        family, standardizer, detector, type, gamma, weight_name,
        bandwidth_h, hset_name, scale_weight_name, method_group,
        method_label, method_id, dgp_type, m, T
      )
    ]
  }), use.names = TRUE, fill = TRUE)
}

summarize_power_files <- function(paths, filter_fun, cv_dt) {
  rbindlist(lapply(paths, function(path) {
    dt <- fread(path, select = selected_cols, showProgress = FALSE)
    dt <- filter_fun(dt)
    if (nrow(dt) == 0L) return(NULL)
    dt[, reject := clean_bool(reject)]
    dt[, T := as.character(T)]
    dt <- merge(
      dt,
      cv_dt[, .(dgp_type, m, T, method_id, empirical_cv)],
      by = c("dgp_type", "m", "T", "method_id"),
      all.x = TRUE,
      sort = FALSE
    )
    dt[, size_adjusted_reject := !is.na(empirical_cv) & statistic > empirical_cv]
    dt[, delay := fifelse(
      reject & !is.na(first_rejection) & first_rejection >= s_star,
      as.numeric(first_rejection - s_star),
      NA_real_
    )]
    dt[
      ,
      .(
        power = mean(reject),
        size_adjusted_power = mean(size_adjusted_reject),
        adl = if (all(is.na(delay))) NA_real_ else mean(delay, na.rm = TRUE),
        n_delay = sum(!is.na(delay)),
        false_alarm_rate = mean(reject & !is.na(first_rejection) & first_rejection < s_star)
      ),
      by = .(
        family, standardizer, detector, type, gamma, weight_name,
        bandwidth_h, hset_name, scale_weight_name, method_group,
        method_label, method_id, scenario, dgp_type, m, T, s_star, delta
      )
    ]
  }), use.names = TRUE, fill = TRUE)
}

aggregate_metric_dt <- function(dt, metric, by_cols) {
  dt[
    ,
    .(
      value = switch(
        metric,
        power = mean(power, na.rm = TRUE),
        sap = mean(size_adjusted_power, na.rm = TRUE),
        adl = if (all(is.na(adl))) NA_real_ else mean(adl, na.rm = TRUE),
        stop("Unknown metric: ", metric)
      )
    ),
    by = by_cols
  ]
}

make_panel_table_lines <- function(data, columns, header, caption, label, formatter, resize = TRUE) {
  tbl <- copy(data)
  for (col in columns) tbl[[col]] <- formatter(tbl[[col]])
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\small",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label)
  )
  if (resize) lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
  lines <- c(
    lines,
    sprintf("\\begin{tabular}{%s}", paste0("l", paste(rep("r", length(columns)), collapse = ""))),
    "\\toprule",
    header,
    "\\midrule"
  )
  panels <- unique(tbl$panel)
  for (i in seq_along(panels)) {
    sub <- tbl[panel == panels[i]]
    lines <- c(lines, sprintf("\\multicolumn{%d}{c}{%s} \\\\", length(columns) + 1L, panels[i]), "\\midrule")
    for (j in seq_len(nrow(sub))) {
      row_vals <- c(sub$paper_label[j], unlist(sub[j, ..columns], use.names = FALSE))
      lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
    }
    if (i < length(panels)) lines <- c(lines, "\\midrule")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  if (resize) lines <- c(lines, "}")
  lines <- c(lines, "\\endgroup", "\\end{table}")
  lines
}

null_overview_table_lines <- function(null_dt, caption, label) {
  tbl <- null_dt[
    ,
    .(
      BB = mean(size[dgp_type == "BB"]),
      fIID = mean(size[dgp_type == "fIID"]),
      fMA1 = mean(size[dgp_type == "fMA1"]),
      `T>=5` = mean(size[T %in% c("5", "10")]),
      MAD = mean(abs(size - 0.05))
    ),
    by = .(panel, paper_label, method_rank)
  ]
  setorder(tbl, method_rank)
  tbl_fmt <- copy(tbl)
  tbl_fmt[, `:=`(
    BB = fmt_pct(BB),
    fIID = fmt_pct(fIID),
    fMA1 = fmt_pct(fMA1),
    `T>=5` = fmt_pct(`T>=5`),
    MAD = fmt_pp(MAD)
  )]
  make_panel_table_lines(
    tbl_fmt,
    columns = c("BB", "fIID", "fMA1", "T>=5", "MAD"),
    header = "Method & BB & fIID & fMA(1) & $T\\ge 5$ & MAD from 5\\% \\\\",
    caption = caption,
    label = label,
    formatter = identity
  )
}

delta_panel_tabular <- function(panel_dt, panel_title, formatter) {
  panel_dt <- copy(panel_dt)
  panel_dt[, delta_label := sprintf("%.3f", delta)]
  wide <- dcast(panel_dt, row_rank + paper_label ~ delta_label, value.var = "value")
  setorder(wide, row_rank)
  cols <- setdiff(names(wide), c("row_rank", "paper_label"))

  lines <- c(
    sprintf("\\textbf{%s}", panel_title),
    "\\smallskip",
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", paste0("l", paste(rep("r", length(cols)), collapse = ""))),
    "\\toprule",
    paste(c("Method", cols), collapse = " & "),
    " \\\\",
    "\\midrule"
  )
  for (j in seq_len(nrow(wide))) {
    row_vals <- c(
      wide$paper_label[j],
      formatter(unlist(wide[j, ..cols], use.names = FALSE))
    )
    lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\par\\smallskip"
  )
  lines
}

delta_table_lines <- function(long_dt, caption, label, value_type = c("pct", "num")) {
  value_type <- match.arg(value_type)
  formatter <- if (identical(value_type, "pct")) fmt_pct else fmt_num
  dgp_titles <- c(
    BB = "Panel A: BB errors",
    fIID = "Panel B: fIID errors",
    fMA1 = "Panel C: fMA(1) errors"
  )

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\tiny",
    "\\setlength{\\tabcolsep}{3.0pt}",
    "\\renewcommand{\\arraystretch}{0.88}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label)
  )

  for (dgp in dgp_levels) {
    sub <- copy(long_dt[dgp_type == dgp])
    if (nrow(sub) == 0L) next
    lines <- c(lines, delta_panel_tabular(sub, dgp_titles[[dgp]], formatter))
  }

  if (tail(lines, 1L) == "\\par\\medskip") {
    lines <- lines[-length(lines)]
  }
  lines <- c(lines, "\\endgroup", "\\end{table}")
  lines
}

aggregate_over_sstar <- function(dt, metric, by_cols) {
  dt[
    ,
    .(
      value = switch(
        metric,
        power = mean(power, na.rm = TRUE),
        sap = mean(size_adjusted_power, na.rm = TRUE),
        adl = if (all(is.na(adl))) NA_real_ else mean(adl, na.rm = TRUE),
        stop("Unknown metric: ", metric)
      )
    ),
    by = by_cols
  ]
}

main_design_table_dt <- function(dt, scenario_value, type_value, metric) {
  sub <- copy(dt[scenario == scenario_value & type == type_value])
  if (metric %in% c("power", "sap", "adl")) {
    sub <- sub[delta > 0]
  }
  sub <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "paper_label", "method_rank", "T", "delta")
  )
  setnames(sub, "method_rank", "row_rank")
  sub[, T := factor(as.character(T), levels = t_levels)]
  sub
}

cvm_design_table_dt <- function(dt, scenario_value, metric) {
  sub <- copy(dt[scenario == scenario_value & type == "CvM"])
  if (metric %in% c("power", "sap", "adl")) {
    sub <- sub[delta > 0]
  }
  sub <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "standardizer", "weight_name", "T", "delta")
  )
  sub[, std_rank := match(standardizer, standardizer_levels)]
  sub[, weight_rank := match(weight_name, weight_levels)]
  sub[, T := factor(as.character(T), levels = t_levels)]
  setorder(sub, dgp_type, weight_rank, std_rank, T, delta)
  sub
}

alt_design_table_dt <- function(dt, scenario_value, metric) {
  sub <- copy(dt[scenario == scenario_value])
  if (metric %in% c("power", "sap", "adl")) {
    sub <- sub[delta > 0]
  }
  sub <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "paper_label", "method_rank", "T", "delta")
  )
  setnames(sub, "method_rank", "row_rank")
  sub[, T := factor(as.character(T), levels = t_levels)]
  sub
}

design_panel_tabular <- function(panel_dt, panel_title, formatter) {
  deltas <- sort(unique(panel_dt$delta))
  Ts <- t_levels[t_levels %in% as.character(unique(panel_dt$T))]
  wide <- dcast(
    copy(panel_dt)[, T := as.character(T)],
    row_rank + paper_label ~ T + delta,
    value.var = "value"
  )
  setorder(wide, row_rank)

  column_keys <- unlist(lapply(Ts, function(Tv) paste0(Tv, "_", deltas)), use.names = FALSE)
  header_top <- c("Method")
  for (Tv in Ts) {
    header_top <- c(header_top, sprintf("\\multicolumn{%d}{c}{$T=%s$}", length(deltas), Tv))
  }

  cmidrules <- character()
  start_col <- 2L
  for (Tv in Ts) {
    end_col <- start_col + length(deltas) - 1L
    cmidrules <- c(cmidrules, sprintf("\\cmidrule(lr){%d-%d}", start_col, end_col))
    start_col <- end_col + 1L
  }

  header_bottom <- c("Method")
  for (Tv in Ts) {
    header_bottom <- c(header_bottom, vapply(deltas, fmt_delta, character(1L)))
  }

  lines <- c(
    sprintf("\\textbf{%s}", panel_title),
    "\\smallskip",
    "\\resizebox{0.97\\textheight}{!}{%",
    sprintf("\\begin{tabular}{%s}", paste0("l", paste(rep("r", length(column_keys)), collapse = ""))),
    "\\toprule",
    paste(header_top, collapse = " & "),
    " \\\\",
    paste(cmidrules, collapse = " "),
    paste(header_bottom, collapse = " & "),
    " \\\\",
    "\\midrule"
  )

  for (j in seq_len(nrow(wide))) {
    row_vals <- c(wide$paper_label[j], formatter(unlist(wide[j, ..column_keys], use.names = FALSE)))
    lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\par\\smallskip"
  )
  lines
}

design_sidewaystable_lines <- function(long_dt, caption, label, value_type = c("pct", "num")) {
  value_type <- match.arg(value_type)
  formatter <- if (identical(value_type, "pct")) fmt_pct else fmt_num
  dgp_titles <- c(
    BB = "Panel A: BB errors",
    fIID = "Panel B: fIID errors",
    fMA1 = "Panel C: fMA(1) errors"
  )

  lines <- c(
    "\\begin{sidewaystable}[p]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\tiny",
    "\\setlength{\\tabcolsep}{2.6pt}",
    "\\renewcommand{\\arraystretch}{0.78}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label)
  )

  for (dgp in dgp_levels) {
    sub <- copy(long_dt[dgp_type == dgp])
    if (nrow(sub) == 0L) next
    lines <- c(lines, design_panel_tabular(sub, dgp_titles[[dgp]], formatter))
  }

  if (length(lines) > 0L && tail(lines, 1L) == "\\par\\smallskip") {
    lines <- lines[-length(lines)]
  }
  lines <- c(lines, "\\endgroup", "\\end{sidewaystable}")
  lines
}

cvm_grouped_rows <- function(wide, value_cols, formatter) {
  lines <- character()
  used_weights <- weight_levels[weight_levels %in% unique(wide$weight_name)]
  for (idx in seq_along(used_weights)) {
    w <- used_weights[idx]
    block <- copy(wide[weight_name == w])
    if (!nrow(block)) next
    setorder(block, std_rank)
    block_label <- unname(weight_table_pretty[[w]])
    for (j in seq_len(nrow(block))) {
      weight_cell <- if (j == 1L) sprintf("\\multirow{%d}{*}{%s}", nrow(block), block_label) else ""
      row_vals <- c(
        weight_cell,
        block$standardizer[j],
        formatter(unlist(block[j, ..value_cols], use.names = FALSE))
      )
      lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
    }
    if (idx < length(used_weights)) lines <- c(lines, "\\addlinespace[2pt]")
  }
  lines
}

design_panel_tabular_cvm <- function(panel_dt, panel_title, formatter) {
  deltas <- sort(unique(panel_dt$delta))
  Ts <- t_levels[t_levels %in% as.character(unique(panel_dt$T))]
  wide <- dcast(
    copy(panel_dt)[, T := as.character(T)],
    weight_rank + weight_name + std_rank + standardizer ~ T + delta,
    value.var = "value"
  )
  setorder(wide, weight_rank, std_rank)

  column_keys <- unlist(lapply(Ts, function(Tv) paste0(Tv, "_", deltas)), use.names = FALSE)
  header_top <- c("Weight", "Method")
  for (Tv in Ts) {
    header_top <- c(header_top, sprintf("\\multicolumn{%d}{c}{$T=%s$}", length(deltas), Tv))
  }

  cmidrules <- character()
  start_col <- 3L
  for (Tv in Ts) {
    end_col <- start_col + length(deltas) - 1L
    cmidrules <- c(cmidrules, sprintf("\\cmidrule(lr){%d-%d}", start_col, end_col))
    start_col <- end_col + 1L
  }

  header_bottom <- c("Weight", "Method")
  for (Tv in Ts) {
    header_bottom <- c(header_bottom, vapply(deltas, fmt_delta, character(1L)))
  }

  lines <- c(
    sprintf("\\textbf{%s}", panel_title),
    "\\smallskip",
    "\\resizebox{0.99\\textheight}{!}{%",
    sprintf("\\begin{tabular}{%s}", paste0("ll", paste(rep("r", length(column_keys)), collapse = ""))),
    "\\toprule",
    paste(header_top, collapse = " & "),
    " \\\\",
    paste(cmidrules, collapse = " "),
    paste(header_bottom, collapse = " & "),
    " \\\\",
    "\\midrule",
    cvm_grouped_rows(wide, column_keys, formatter),
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\par\\smallskip"
  )
  lines
}

design_sidewaystable_lines_cvm <- function(long_dt, caption, label, value_type = c("pct", "num")) {
  value_type <- match.arg(value_type)
  formatter <- if (identical(value_type, "pct")) fmt_pct else fmt_num
  dgp_titles <- c(
    BB = "Panel A: BB errors",
    fIID = "Panel B: fIID errors",
    fMA1 = "Panel C: fMA(1) errors"
  )

  lines <- c(
    "\\begin{sidewaystable}[p]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\tiny",
    "\\setlength{\\tabcolsep}{2.2pt}",
    "\\renewcommand{\\arraystretch}{0.76}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label)
  )

  for (dgp in dgp_levels) {
    sub <- copy(long_dt[dgp_type == dgp])
    if (nrow(sub) == 0L) next
    lines <- c(lines, design_panel_tabular_cvm(sub, dgp_titles[[dgp]], formatter))
  }

  if (length(lines) > 0L && tail(lines, 1L) == "\\par\\smallskip") {
    lines <- lines[-length(lines)]
  }
  lines <- c(lines, "\\endgroup", "\\end{sidewaystable}")
  lines
}

null_detail_table_core_lines <- function(sub, caption, label, m_value = 500L) {
  sub <- copy(sub[m == m_value])
  sub[, T := factor(as.character(T), levels = t_levels)]
  wide <- dcast(
    sub,
    paper_label + method_rank ~ dgp_type + T,
    value.var = "size"
  )
  setorder(wide, method_rank)

  dgp_used <- dgp_levels[dgp_levels %in% unique(sub$dgp_type)]
  cols <- unlist(lapply(dgp_used, function(dg) paste(dg, t_levels, sep = "_")), use.names = FALSE)

  top <- c("Method")
  for (dg in dgp_used) {
    top <- c(top, sprintf("\\multicolumn{%d}{c}{%s}", length(t_levels), dgp_pretty[[dg]]))
  }
  cmidrules <- character()
  start_col <- 2L
  for (dg in dgp_used) {
    end_col <- start_col + length(t_levels) - 1L
    cmidrules <- c(cmidrules, sprintf("\\cmidrule(lr){%d-%d}", start_col, end_col))
    start_col <- end_col + 1L
  }
  bottom <- c("Method", rep(paste0("$T=", t_levels, "$"), times = length(dgp_used)))

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\scriptsize",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", paste0("l", paste(rep("r", length(cols)), collapse = ""))),
    "\\toprule",
    paste(top, collapse = " & "),
    " \\\\",
    paste(cmidrules, collapse = " "),
    paste(bottom, collapse = " & "),
    " \\\\",
    "\\midrule"
  )

  for (j in seq_len(nrow(wide))) {
    row_vals <- c(wide$paper_label[j], fmt_pct(unlist(wide[j, ..cols], use.names = FALSE)))
    lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\endgroup",
    "\\end{table}"
  )
  lines
}

null_detail_table_lines <- function(null_dt, type_value, caption, label, m_value = 500L) {
  null_detail_table_core_lines(
    sub = null_dt[type == type_value],
    caption = caption,
    label = label,
    m_value = m_value
  )
}

null_detail_table_lines_cvm <- function(null_dt, caption, label, m_value = 500L) {
  sub <- copy(null_dt[type == "CvM" & m == m_value])
  sub[, std_rank := match(standardizer, standardizer_levels)]
  sub[, weight_rank := match(weight_name, weight_levels)]
  sub[, T := factor(as.character(T), levels = t_levels)]
  wide <- dcast(
    sub,
    weight_rank + weight_name + std_rank + standardizer ~ dgp_type + T,
    value.var = "size"
  )
  setorder(wide, weight_rank, std_rank)

  dgp_used <- dgp_levels[dgp_levels %in% unique(sub$dgp_type)]
  cols <- unlist(lapply(dgp_used, function(dg) paste(dg, t_levels, sep = "_")), use.names = FALSE)

  header_top <- c("Weight", "Method")
  for (dg in dgp_used) {
    header_top <- c(header_top, sprintf("\\multicolumn{%d}{c}{%s}", length(t_levels), dgp_pretty[[dg]]))
  }
  cmidrules <- character()
  start_col <- 3L
  for (dg in dgp_used) {
    end_col <- start_col + length(t_levels) - 1L
    cmidrules <- c(cmidrules, sprintf("\\cmidrule(lr){%d-%d}", start_col, end_col))
    start_col <- end_col + 1L
  }
  header_bottom <- c("Weight", "Method")
  for (dg in dgp_used) header_bottom <- c(header_bottom, paste0("$T=", t_levels, "$"))

  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\singlespacing",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2.8pt}",
    "\\renewcommand{\\arraystretch}{0.84}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", paste0("ll", paste(rep("r", length(cols)), collapse = ""))),
    "\\toprule",
    paste(header_top, collapse = " & "),
    " \\\\",
    paste(cmidrules, collapse = " "),
    paste(header_bottom, collapse = " & "),
    " \\\\",
    "\\midrule",
    cvm_grouped_rows(wide, cols, fmt_pct),
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\endgroup",
    "\\end{table}"
  )
}

subfigure_lines <- function(path, caption, width) {
  c(
    sprintf("\\begin{subfigure}[t]{%s}", width),
    "\\centering",
    sprintf("\\includegraphics[width=\\linewidth]{%s}", path),
    sprintf("\\caption{%s}", caption),
    "\\end{subfigure}"
  )
}

subfigure_grid_lines <- function(paths, captions, widths, per_row = 2L) {
  out <- character()
  n <- length(paths)
  for (i in seq_len(n)) {
    out <- c(out, subfigure_lines(paths[i], captions[i], widths[i]))
    if (i < n) {
      if (i %% per_row == 0L) {
        out <- c(out, "\\par\\medskip")
      } else {
        out <- c(out, "\\hfill")
      }
    }
  }
  out
}

figure_with_subfigures <- function(paths, captions, widths, figure_caption, label, per_row = 2L) {
  c(
    "\\begin{figure}[!htbp]",
    "\\centering",
    subfigure_grid_lines(paths, captions, widths, per_row = per_row),
    sprintf("\\caption{%s}", figure_caption),
    sprintf("\\label{%s}", label),
    "\\end{figure}"
  )
}

three_dgp_figure_lines <- function(paths, figure_caption, label) {
  figure_with_subfigures(
    paths = paths,
    captions = vapply(dgp_levels, dgp_subcaption, character(1L)),
    widths = c("0.84\\textwidth", "0.84\\textwidth", "0.84\\textwidth"),
    figure_caption = figure_caption,
    label = label,
    per_row = 1L
  )
}

metric_curve_caption <- function(metric) {
  switch(
    metric,
    power = "Raw rejection curves",
    sap = "Size-adjusted rejection curves",
    adl = "Average detection lag curves"
  )
}

ks_style_caption <- paste(
  "Within every panel, RSMS with $\\gamma=0$ is the black solid line with filled circles,",
  "RSMS with $\\gamma=0.15$ is the black long-dashed line with open circles,",
  "SSMS with $\\gamma=0$ is the dark-gray dashed line with filled squares,",
  "SSMS with $\\gamma=0.15$ is the dark-gray dot-dashed line with open squares,",
  "HAC with $\\gamma=0$ is the light-gray dotted line with filled triangles,",
  "and HAC with $\\gamma=0.15$ is the light-gray two-dash line with open triangles."
)

std_style_caption <- paste(
  "Within every panel, RSMS is the black solid line with circles,",
  "SSMS is the dark-gray dashed line with squares,",
  "and HAC is the light-gray dotted line with triangles."
)

cvm_grouped_table_note <- paste(
  "Within each DGP panel, rows are grouped by the Uniform, Early, Mid, and Late weights,",
  "and each weight block reports RSMS, SSMS, and HAC."
)

replace_between_markers <- function(lines, start_marker, end_marker, new_lines) {
  start_idx <- grep(start_marker, lines, fixed = TRUE)
  end_idx <- grep(end_marker, lines, fixed = TRUE)
  if (length(start_idx) != 1L || length(end_idx) != 1L || end_idx <= start_idx) {
    stop("Could not locate marker pair: ", start_marker, " / ", end_marker)
  }
  c(
    lines[seq_len(start_idx)],
    new_lines,
    lines[end_idx:length(lines)]
  )
}

plot_theme <- function(base_size = 8.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.28),
      panel.spacing.x = grid::unit(0.55, "lines"),
      strip.background = element_rect(fill = "grey94", colour = "grey65"),
      strip.text = element_text(size = base_size - 1, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

ks_series_levels <- c(
  "RSMS, gamma=0",
  "RSMS, gamma=0.15",
  "SSMS, gamma=0",
  "SSMS, gamma=0.15",
  "HAC, gamma=0",
  "HAC, gamma=0.15"
)
ks_colors <- c(
  "RSMS, gamma=0" = "black",
  "RSMS, gamma=0.15" = "black",
  "SSMS, gamma=0" = "grey35",
  "SSMS, gamma=0.15" = "grey35",
  "HAC, gamma=0" = "grey60",
  "HAC, gamma=0.15" = "grey60"
)
ks_linetypes <- c(
  "RSMS, gamma=0" = "solid",
  "RSMS, gamma=0.15" = "longdash",
  "SSMS, gamma=0" = "dashed",
  "SSMS, gamma=0.15" = "dotdash",
  "HAC, gamma=0" = "dotted",
  "HAC, gamma=0.15" = "twodash"
)
ks_shapes <- c(
  "RSMS, gamma=0" = 16,
  "RSMS, gamma=0.15" = 1,
  "SSMS, gamma=0" = 15,
  "SSMS, gamma=0.15" = 0,
  "HAC, gamma=0" = 17,
  "HAC, gamma=0.15" = 2
)

std_colors <- c(RSMS = "black", SSMS = "grey35", HAC = "grey60")
std_linetypes <- c(RSMS = "solid", SSMS = "dashed", HAC = "dotted")
std_shapes <- c(RSMS = 16, SSMS = 15, HAC = 17)

metric_label <- function(metric) {
  switch(
    metric,
    power = "Raw rejection probability",
    sap = "Size-adjusted rejection probability",
    adl = "Average detection lag"
  )
}

save_plot <- function(p, filename, width, height) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

make_ks_plot_data <- function(dt, scenario_value, metric) {
  sub <- copy(dt[scenario == scenario_value & type == "KS" & delta > 0])
  out <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "T", "delta", "standardizer", "gamma")
  )
  out[, series := factor(
    sprintf("%s, gamma=%s", standardizer, vapply(gamma, gamma_label, character(1L))),
    levels = ks_series_levels
  )]
  out[, dgp_type := factor(dgp_type, levels = dgp_levels, labels = unname(dgp_pretty[dgp_levels]))]
  out[, T := factor(as.character(T), levels = t_levels, labels = paste0("T=", t_levels))]
  out
}

make_cvm_plot_data <- function(dt, scenario_value, weight_value, metric) {
  sub <- copy(dt[scenario == scenario_value & type == "CvM" & weight_name == weight_value & delta > 0])
  out <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "T", "delta", "standardizer")
  )
  out[, standardizer := factor(standardizer, levels = standardizer_levels)]
  out[, dgp_type := factor(dgp_type, levels = dgp_levels, labels = unname(dgp_pretty[dgp_levels]))]
  out[, T := factor(as.character(T), levels = t_levels, labels = paste0("T=", t_levels))]
  out
}

make_alt_plot_data <- function(dt, scenario_value, tuning_value, metric) {
  sub <- copy(dt[scenario == scenario_value & tuning_label_plot == tuning_value & delta > 0])
  out <- aggregate_over_sstar(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "T", "delta", "standardizer")
  )
  out[, standardizer := factor(standardizer, levels = standardizer_levels)]
  out[, dgp_type := factor(dgp_type, levels = dgp_levels, labels = unname(dgp_pretty[dgp_levels]))]
  out[, T := factor(as.character(T), levels = t_levels, labels = paste0("T=", t_levels))]
  out
}

save_ks_scenario_plot <- function(dt, scenario_value, metric, filename) {
  plot_dt <- make_ks_plot_data(dt, scenario_value, metric)
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = series, colour = series, linetype = series, shape = series)
  ) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.15, stroke = 0.25) +
    facet_grid(dgp_type ~ T, scales = "free_x") +
    scale_color_manual(values = ks_colors, breaks = ks_series_levels) +
    scale_linetype_manual(values = ks_linetypes, breaks = ks_series_levels) +
    scale_shape_manual(values = ks_shapes, breaks = ks_series_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    plot_theme(base_size = 8.0)
  save_plot(p, filename, width = 8.6, height = 5.0)
}

metric_axis_limits <- function(metric) {
  switch(
    metric,
    power = coord_cartesian(ylim = c(0, 1.02)),
    sap = coord_cartesian(ylim = c(0, 1.02)),
    adl = NULL
  )
}

dgp_slug <- function(x) {
  out <- c(BB = "bb", fIID = "fiid", fMA1 = "fma1")
  unname(out[x])
}

dgp_subcaption <- function(dgp_value) {
  sprintf(
    "%s errors; the four internal columns, from left to right, correspond to $T=1,2,5,10$.",
    dgp_pretty[[dgp_value]]
  )
}

save_ks_dgp_plot <- function(dt, scenario_value, dgp_value, metric, filename) {
  plot_dt <- make_ks_plot_data(dt, scenario_value, metric)[as.character(dgp_type) == dgp_pretty[[dgp_value]]]
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = series, colour = series, linetype = series, shape = series)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.2, stroke = 0.25) +
    facet_wrap(~T, nrow = 1, scales = "free_x") +
    scale_color_manual(values = ks_colors, breaks = ks_series_levels) +
    scale_linetype_manual(values = ks_linetypes, breaks = ks_series_levels) +
    scale_shape_manual(values = ks_shapes, breaks = ks_series_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    plot_theme(base_size = 8.2) +
    theme(legend.position = "none")
  if (!is.null(metric_axis_limits(metric))) p <- p + metric_axis_limits(metric)
  save_plot(p, filename, width = 7.2, height = 2.45)
}

save_cvm_weight_plot <- function(dt, scenario_value, weight_value, metric, filename) {
  plot_dt <- make_cvm_plot_data(dt, scenario_value, weight_value, metric)
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = standardizer, colour = standardizer, linetype = standardizer, shape = standardizer)
  ) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.15, stroke = 0.25) +
    facet_grid(dgp_type ~ T, scales = "free_x") +
    scale_color_manual(values = std_colors, breaks = standardizer_levels) +
    scale_linetype_manual(values = std_linetypes, breaks = standardizer_levels) +
    scale_shape_manual(values = std_shapes, breaks = standardizer_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    plot_theme(base_size = 8.0)
  save_plot(p, filename, width = 8.6, height = 5.0)
}

save_cvm_dgp_plot <- function(dt, scenario_value, weight_value, dgp_value, metric, filename) {
  plot_dt <- make_cvm_plot_data(dt, scenario_value, weight_value, metric)[
    as.character(dgp_type) == dgp_pretty[[dgp_value]]
  ]
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = standardizer, colour = standardizer, linetype = standardizer, shape = standardizer)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.2, stroke = 0.25) +
    facet_wrap(~T, nrow = 1, scales = "free_x") +
    scale_color_manual(values = std_colors, breaks = standardizer_levels) +
    scale_linetype_manual(values = std_linetypes, breaks = standardizer_levels) +
    scale_shape_manual(values = std_shapes, breaks = standardizer_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    plot_theme(base_size = 8.2) +
    theme(legend.position = "none")
  if (!is.null(metric_axis_limits(metric))) p <- p + metric_axis_limits(metric)
  save_plot(p, filename, width = 7.2, height = 2.45)
}

save_alt_tuning_plot <- function(dt, scenario_value, tuning_value, metric, filename) {
  plot_dt <- make_alt_plot_data(dt, scenario_value, tuning_value, metric)
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = standardizer, colour = standardizer, linetype = standardizer, shape = standardizer)
  ) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.15, stroke = 0.25) +
    facet_grid(dgp_type ~ T, scales = "free_x") +
    scale_color_manual(values = std_colors, breaks = standardizer_levels) +
    scale_linetype_manual(values = std_linetypes, breaks = standardizer_levels) +
    scale_shape_manual(values = std_shapes, breaks = standardizer_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    plot_theme(base_size = 8.0)
  save_plot(p, filename, width = 8.6, height = 5.0)
}

save_alt_dgp_plot <- function(dt, scenario_value, tuning_value, dgp_value, metric, filename) {
  plot_dt <- make_alt_plot_data(dt, scenario_value, tuning_value, metric)[
    as.character(dgp_type) == dgp_pretty[[dgp_value]]
  ]
  p <- ggplot(
    plot_dt,
    aes(x = delta, y = value, group = standardizer, colour = standardizer, linetype = standardizer, shape = standardizer)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.2, stroke = 0.25) +
    facet_wrap(~T, nrow = 1, scales = "free_x") +
    scale_color_manual(values = std_colors, breaks = standardizer_levels) +
    scale_linetype_manual(values = std_linetypes, breaks = standardizer_levels) +
    scale_shape_manual(values = std_shapes, breaks = standardizer_levels) +
    labs(x = expression(Delta), y = metric_label(metric)) +
    plot_theme(base_size = 8.2) +
    theme(legend.position = "none")
  if (!is.null(metric_axis_limits(metric))) p <- p + metric_axis_limits(metric)
  save_plot(p, filename, width = 7.2, height = 2.45)
}

main_long_table_dt <- function(dt, scenario_value, type_value, metric) {
  sub <- copy(dt[scenario == scenario_value & type == type_value])
  if (identical(metric, "adl")) sub <- sub[delta > 0]
  else sub <- sub[delta > 0]
  out <- aggregate_metric_dt(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "paper_label", "method_rank", "delta")
  )
  setnames(out, "method_rank", "row_rank")
  out
}

alt_long_table_dt <- function(dt, scenario_value, metric) {
  sub <- copy(dt[scenario == scenario_value])
  if (identical(metric, "adl")) sub <- sub[delta > 0]
  else sub <- sub[delta > 0]
  out <- aggregate_metric_dt(
    sub,
    metric = metric,
    by_cols = c("dgp_type", "paper_label", "method_rank", "delta")
  )
  setnames(out, "method_rank", "row_rank")
  out
}

metric_table_caption <- function(metric, subject_text, scenario_value) {
  metric_text <- switch(
    metric,
    power = "Raw rejection rates (in \\%)",
    sap = "Size-adjusted rejection rates (in \\%)",
    adl = "Average detection lag"
  )
  tail_text <- switch(
    metric,
    power = "Each entry averages over the available monitoring horizons and break locations while keeping the DGP-specific positive break magnitudes $\\Delta$ fixed.",
    sap = "Each entry averages over the available monitoring horizons and break locations while recalibrating by the empirical null critical value from the matching design.",
    adl = "Each entry averages over the available monitoring horizons and break locations at a fixed positive break magnitude, and smaller values indicate faster post-break detection."
  )
  sprintf("%s for %s under the %s setting. %s", metric_text, subject_text, scenario_pretty_adj[[scenario_value]], tail_text)
}

metric_design_table_caption <- function(metric, subject_text, scenario_value) {
  metric_text <- switch(
    metric,
    power = "Raw rejection rates (in \\%)",
    sap = "Size-adjusted rejection rates (in \\%)",
    adl = "Average detection lag"
  )
  tail_text <- switch(
    metric,
    power = "Panels separate the three DGPs, columns are grouped by monitoring horizon $T\\in\\{1,2,5,10\\}$ and then by the positive break sizes $\\Delta$, and entries average over the two break locations $s^{\\star}\\in\\{50,200\\}$.",
    sap = "Panels separate the three DGPs, columns are grouped by monitoring horizon $T\\in\\{1,2,5,10\\}$ and then by the positive break sizes $\\Delta$, and the empirical size adjustment is computed design by design before averaging over $s^{\\star}\\in\\{50,200\\}$.",
    adl = "Panels separate the three DGPs, columns are grouped by monitoring horizon $T\\in\\{1,2,5,10\\}$ and then by the positive break sizes $\\Delta$, and smaller entries indicate faster post-break detection after averaging over $s^{\\star}\\in\\{50,200\\}$."
  )
  sprintf("%s for %s under the %s setting. %s", metric_text, subject_text, scenario_pretty_adj[[scenario_value]], tail_text)
}

fmt_pct_text <- function(x) {
  if (!length(x) || !is.finite(x[1L])) return("NA")
  sprintf("%.1f\\%%", 100 * x[1L])
}

fmt_num_text <- function(x) {
  if (!length(x) || !is.finite(x[1L])) return("NA")
  sprintf("%.1f", x[1L])
}

pct_range_text <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf("%s--%s", fmt_pct_text(min(x)), fmt_pct_text(max(x)))
}

num_range_text <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf("%s--%s", fmt_num_text(min(x)), fmt_num_text(max(x)))
}

format_t_set <- function(x) {
  vals <- sort(unique(as.numeric(as.character(x))))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return("$T\\in\\{\\}$")
  sprintf(
    "$T\\in\\{%s\\}$",
    paste(vapply(vals, function(v) formatC(v, format = "f", digits = 0), character(1L)), collapse = ",")
  )
}

oxford_join <- function(x) {
  x <- x[nzchar(x)]
  if (!length(x)) return("")
  if (length(x) == 1L) return(x)
  if (length(x) == 2L) return(sprintf("%s and %s", x[1L], x[2L]))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}

metric_ordered <- function(dt, metric) {
  if (!nrow(dt)) return(dt)
  if (identical(metric, "adl")) {
    dt[order(value)]
  } else {
    dt[order(-value)]
  }
}

aggregate_metric_nonnull <- function(dt, metric, by_cols) {
  aggregate_metric_dt(dt[delta > 0], metric, by_cols)
}

top_metric_row <- function(dt, metric, by_cols) {
  agg <- aggregate_metric_nonnull(dt, metric, by_cols)
  agg <- agg[is.finite(value)]
  if (!nrow(agg)) return(NULL)
  metric_ordered(agg, metric)[1L]
}

bottom_metric_row <- function(dt, metric, by_cols) {
  agg <- aggregate_metric_nonnull(dt, metric, by_cols)
  agg <- agg[is.finite(value)]
  if (!nrow(agg)) return(NULL)
  if (identical(metric, "adl")) {
    agg[order(-value)][1L]
  } else {
    agg[order(value)][1L]
  }
}

winner_by_t <- function(dt, metric, by_cols) {
  agg <- aggregate_metric_nonnull(dt, metric, c("T", by_cols))
  agg <- agg[is.finite(value)]
  if (!nrow(agg)) return(agg)
  if (identical(metric, "adl")) {
    agg[, .SD[which.min(value)][1L], by = T]
  } else {
    agg[, .SD[which.max(value)][1L], by = T]
  }
}

dominant_winner_phrase <- function(winner_dt, col, pretty_map = NULL, verb = NULL) {
  vals <- winner_dt[[col]]
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return("no stable winner")
  uniq <- unique(vals)
  if (length(uniq) == 1L) {
    label <- if (is.null(pretty_map)) uniq[1L] else unname(pretty_map[[uniq[1L]]])
    if (is.null(verb)) {
      return(sprintf("%s at every horizon $T=1,2,5,10$", label))
    }
    return(sprintf("%s %s at every horizon $T=1,2,5,10$", label, verb))
  }
  counts <- sort(table(vals), decreasing = TRUE)
  lead <- names(counts)[1L]
  lead_label <- if (is.null(pretty_map)) lead else unname(pretty_map[[lead]])
  lead_t <- winner_dt$T[vals == lead]
  if (is.null(verb)) {
    return(sprintf("%s on most horizons (specifically %s)", lead_label, format_t_set(lead_t)))
  }
  sprintf("%s %s on most horizons (specifically %s)", lead_label, verb, format_t_set(lead_t))
}

ordered_group_labels <- function(dt, metric, group_col, pretty_map = NULL) {
  agg <- aggregate_metric_nonnull(dt, metric, by_cols = group_col)
  agg <- agg[is.finite(value)]
  if (!nrow(agg)) return(character())
  vals <- metric_ordered(agg, metric)[[group_col]]
  if (!is.null(pretty_map)) vals <- unname(pretty_map[vals])
  vals
}

metric_path_t <- function(dt, metric) {
  out <- aggregate_metric_nonnull(dt, metric, by_cols = "T")
  out[, T := as.character(T)]
  out
}

metric_value_t <- function(path_dt, T_value) {
  path_dt[as.character(T) == as.character(T_value), value][1L]
}

progress_sentence <- function(path_dt, metric, subject_label, scenario_value = NULL) {
  v1 <- metric_value_t(path_dt, "1")
  v2 <- metric_value_t(path_dt, "2")
  v5 <- metric_value_t(path_dt, "5")
  v10 <- metric_value_t(path_dt, "10")
  if (identical(metric, "power")) {
    base <- sprintf(
      "For %s, the mean raw rejection rate rises from %s at $T=1$ to %s at $T=10$.",
      subject_label,
      fmt_pct_text(v1),
      fmt_pct_text(v10)
    )
  } else if (identical(metric, "sap")) {
    base <- sprintf(
      "For %s, the mean size-adjusted rejection rate rises from %s at $T=1$ to %s at $T=10$.",
      subject_label,
      fmt_pct_text(v1),
      fmt_pct_text(v10)
    )
  } else {
    base <- sprintf(
      "For %s, the mean ADL increases from %s at $T=1$ to %s at $T=10$.",
      subject_label,
      fmt_num_text(v1),
      fmt_num_text(v10)
    )
  }
  tail <- ""
  if (identical(metric, "adl")) {
    tail <- "ADL therefore rises in absolute units as the monitoring horizon lengthens, because stopping times are recorded on a longer scale."
  } else if (all(is.finite(c(v1, v2, v5, v10)))) {
    tail <- if ((v2 - v1) >= (v10 - v5)) {
      "Most of the gain materializes between $T=1$ and $T=2$, with smaller increments thereafter."
    } else {
      "The longer horizons continue to deliver visible gains even after $T=5$."
    }
  }
  if (!is.null(scenario_value) && scenario_value == "smooth_change" && metric %in% c("power", "sap")) {
    tail <- paste(
      tail,
      "This horizon effect is especially pronounced for smooth change, where the drift needs time to accumulate."
    )
  }
  paste(base, tail)
}

ks_gamma_sentence <- function(dt, metric) {
  gamma_rank <- top_metric_row(dt, metric, by_cols = "gamma")
  if (is.null(gamma_rank) || !is.finite(gamma_rank$gamma)) return("")
  if (identical(metric, "adl")) {
    if (abs(gamma_rank$gamma - 0.15) < 1e-12) {
      "Within the self-normalized KS class, the $\\gamma=0.15$ boundary is the faster ADL choice."
    } else {
      "Within the self-normalized KS class, the $\\gamma=0$ boundary is the faster ADL choice."
    }
  } else {
    if (abs(gamma_rank$gamma) < 1e-12) {
      "Within the self-normalized KS class, the $\\gamma=0$ boundary is the more powerful choice."
    } else {
      "Within the self-normalized KS class, the $\\gamma=0.15$ boundary is the more powerful choice."
    }
  }
}

rsms_vs_ssms_gap_sentence <- function(dt, metric, best_row) {
  if (is.null(best_row) || !identical(best_row$standardizer, "RSMS") || !is.finite(best_row$gamma)) return("")
  rsms_dt <- dt[standardizer == "RSMS" & abs(gamma - best_row$gamma) < 1e-12]
  ssms_dt <- dt[standardizer == "SSMS" & abs(gamma - best_row$gamma) < 1e-12]
  if (!nrow(rsms_dt) || !nrow(ssms_dt)) return("")
  rsms_path <- metric_path_t(rsms_dt, metric)
  ssms_path <- metric_path_t(ssms_dt, metric)
  r1 <- metric_value_t(rsms_path, "1")
  r10 <- metric_value_t(rsms_path, "10")
  s1 <- metric_value_t(ssms_path, "1")
  s10 <- metric_value_t(ssms_path, "10")
  if (identical(metric, "adl")) {
    if (!all(is.finite(c(r1, r10, s1, s10))) || !(r1 < s1 && r10 < s10)) return("")
    sprintf(
      "Relative to the matching SSMS boundary, RSMS reduces the mean delay by about %s and %s monitoring units at $T=1$ and $T=10$, respectively.",
      fmt_num_text(s1 - r1),
      fmt_num_text(s10 - r10)
    )
  } else {
    if (!all(is.finite(c(r1, r10, s1, s10))) || !(r1 > s1 && r10 > s10)) return("")
    sprintf(
      "Relative to the matching SSMS boundary, RSMS raises the mean rejection rate by about %.1f and %.1f percentage points at $T=1$ and $T=10$, respectively.",
      100 * (r1 - s1),
      100 * (r10 - s10)
    )
  }
}

main_null_ks_comment <- function(null_dt) {
  sub <- null_dt[type == "KS"]
  by_std_t <- sub[, .(size = mean(size, na.rm = TRUE)), by = .(standardizer, T)]
  best_t <- by_std_t[, .SD[which.min(abs(size - 0.05))], by = T]
  best_phrase <- if (length(unique(best_t$standardizer)) == 1L) {
    sprintf("%s is closest to the nominal 5\\%% size at every horizon $T=1,2,5,10$.", unique(best_t$standardizer))
  } else {
    sprintf("%s is closest to nominal on most horizons.", names(sort(table(best_t$standardizer), decreasing = TRUE))[1L])
  }
  rsms_t1 <- by_std_t[standardizer == "RSMS" & as.character(T) == "1", size][1L]
  rsms_t10 <- by_std_t[standardizer == "RSMS" & as.character(T) == "10", size][1L]
  paste(
    "The horizon-specific results in Table~\\ref{tab:sim-null-ks} refine the size comparison.",
    best_phrase,
    sprintf(
      "Averaging over DGPs and the two $\\gamma$ values, SSMS KS stays in the %s band, RSMS ranges from about %s at $T=1$ to %s at $T=10$, and HAC remains the most liberal benchmark at roughly %s.",
      pct_range_text(by_std_t[standardizer == "SSMS", size]),
      fmt_pct_text(rsms_t1),
      fmt_pct_text(rsms_t10),
      pct_range_text(by_std_t[standardizer == "HAC", size])
    ),
    "The subsequent raw-power comparisons should therefore be interpreted jointly with this calibration profile: HAC attains higher rejection rates partly through greater size distortion, whereas RSMS delivers the stronger power-size compromise within the self-normalized class."
  )
}

main_null_cvm_comment <- function(null_dt) {
  sub <- null_dt[type == "CvM"]
  by_std_t <- sub[, .(size = mean(size, na.rm = TRUE)), by = .(standardizer, T)]
  best_t <- by_std_t[, .SD[which.min(abs(size - 0.05))], by = T]
  best_phrase <- if (length(unique(best_t$standardizer)) == 1L) {
    sprintf("%s is again closest to the nominal 5\\%% size at every horizon $T=1,2,5,10$.", unique(best_t$standardizer))
  } else {
    sprintf("%s is closest to nominal on most horizons.", names(sort(table(best_t$standardizer), decreasing = TRUE))[1L])
  }
  ssms_weight_best <- sub[
    standardizer == "SSMS",
    .(mad = mean(abs(size - 0.05), na.rm = TRUE)),
    by = .(weight_name, T)
  ][, .SD[which.min(mad)], by = T]
  paste(
    "Table~\\ref{tab:sim-null-cvm} yields the same conclusion for weighted-CvM monitoring.",
    best_phrase,
    sprintf(
      "Across the four weight families, SSMS weighted-CvM sizes stay in the %s band, RSMS improves from about %s at $T=1$ to %s at $T=10$, and HAC remains materially oversized at roughly %s.",
      pct_range_text(by_std_t[standardizer == "SSMS", size]),
      fmt_pct_text(by_std_t[standardizer == "RSMS" & as.character(T) == "1", size][1L]),
      fmt_pct_text(by_std_t[standardizer == "RSMS" & as.character(T) == "10", size][1L]),
      pct_range_text(by_std_t[standardizer == "HAC", size])
    ),
    sprintf(
      "Within SSMS, %s.",
      dominant_winner_phrase(ssms_weight_best, "weight_name", weight_short_pretty, verb = "is closest to nominal")
    )
  )
}

make_ks_metric_comment <- function(power_dt, scenario_value, metric) {
  sub <- power_dt[scenario == scenario_value & type == "KS" & delta > 0]
  self_sub <- sub[standardizer %in% c("RSMS", "SSMS")]
  overall_best <- top_metric_row(sub, metric, c("paper_label", "standardizer", "gamma"))
  self_best <- top_metric_row(self_sub, metric, c("paper_label", "standardizer", "gamma"))
  if (is.null(overall_best) || is.null(self_best)) return("")
  self_phrase <- if (identical(metric, "adl")) {
    dominant_winner_phrase(winner_by_t(self_sub, metric, c("paper_label")), "paper_label", verb = "has the smallest delay")
  } else {
    dominant_winner_phrase(winner_by_t(self_sub, metric, c("paper_label")), "paper_label", verb = "leads")
  }
  gap_sentence <- rsms_vs_ssms_gap_sentence(self_sub, metric, self_best)
  gamma_sentence <- ks_gamma_sentence(self_sub, metric)
  progress <- progress_sentence(
    metric_path_t(self_sub[paper_label == self_best$paper_label], metric),
    metric,
    self_best$paper_label,
    scenario_value
  )
  if (identical(metric, "power")) {
    paste(
      "The KS raw-power tables and subfigures indicate a stable ranking across monitoring horizons.",
      sprintf(
        "On average, %s attains the largest raw rejection rates, and within the self-normalized class %s.",
        overall_best$paper_label,
        self_phrase
      ),
      gap_sentence,
      gamma_sentence,
      progress
    )
  } else if (identical(metric, "sap")) {
    paste(
      "The size-adjusted KS tables and subfigures provide the corresponding calibration-adjusted comparison.",
      sprintf(
        "%s remains the overall SAP leader, while within the self-normalized class %s once size control is taken into account.",
        overall_best$paper_label,
        self_phrase
      ),
      gap_sentence,
      gamma_sentence,
      progress
    )
  } else {
    paste(
      "The KS ADL tables and subfigures yield the same RSMS-versus-SSMS ordering across horizons.",
      sprintf("Among the self-normalized rules, %s.", self_phrase),
      gap_sentence,
      gamma_sentence,
      progress
    )
  }
}

make_cvm_metric_comment <- function(power_dt, scenario_value, metric) {
  sub <- power_dt[scenario == scenario_value & type == "CvM" & delta > 0]
  rsms_sub <- sub[standardizer == "RSMS"]
  if (!nrow(sub) || !nrow(rsms_sub)) return("")
  rsms_best <- top_metric_row(rsms_sub, metric, c("paper_label", "weight_name"))
  if (is.null(rsms_best)) return("")
  rsms_weight_phrase <- if (identical(metric, "adl")) {
    dominant_winner_phrase(winner_by_t(rsms_sub, metric, c("weight_name")), "weight_name", weight_short_pretty, verb = "has the smallest delay")
  } else {
    dominant_winner_phrase(winner_by_t(rsms_sub, metric, c("weight_name")), "weight_name", weight_short_pretty, verb = "leads")
  }
  order_text <- oxford_join(ordered_group_labels(rsms_sub, metric, "weight_name", weight_short_pretty))
  progress <- progress_sentence(
    metric_path_t(rsms_sub[weight_name == rsms_best$weight_name], metric),
    metric,
    rsms_best$paper_label,
    scenario_value
  )
  if (identical(metric, "adl")) {
    paste(
      "The weighted-CvM ADL tables and subfigures reveal a regular trade-off across weight functions.",
      sprintf("Within RSMS, %s.", rsms_weight_phrase),
      sprintf("Averaging over DGPs and horizons, the RSMS weight ranking is %s.", order_text),
      progress,
      "This pattern highlights the central power-delay trade-off within the weighted-CvM class: earlier weighting reacts more quickly, whereas later weighting waits for a larger accumulation of post-change evidence."
    )
  } else {
    paste(
      sprintf(
        "The weighted-CvM %s tables and subfigures display a similarly regular ordering across weight functions.",
        if (identical(metric, "power")) "raw-power" else "size-adjusted-power"
      ),
      sprintf("Within RSMS, %s.", rsms_weight_phrase),
      sprintf("Averaging over DGPs and horizons, the RSMS weight ranking is %s.", order_text),
      progress,
      "The corresponding SSMS results remain at lower rejection levels, so RSMS emerges as the stronger self-normalized weighted-CvM specification in these designs."
    )
  }
}

alt_null_comment <- function(null_dt) {
  paste(
    "Table~\\ref{tab:app-alt-null-detail} shows that calibration in the alternative detector family is much less uniform than in the main KS/CvM section.",
    sprintf(
      "Within the retained appendix comparators, Page-CUSUM with SSMS stays closest to nominal at roughly %s, whereas weighted CUSUM with HAC reaches %s and multiscale MOSUM with HAC reaches %s.",
      pct_range_text(null_dt[detector == "PageCUSUM" & standardizer == "SSMS", size]),
      pct_range_text(null_dt[detector == "WeightedCUSUM" & standardizer == "HAC", size]),
      pct_range_text(null_dt[detector == "MultiscaleMOSUM" & standardizer == "HAC", size])
    ),
    "SSMS is again the safest standardizer, while HAC and some RSMS combinations are substantially oversized. The current public replication bundle retains the Page-CUSUM, weighted-CUSUM, and multiscale-MOSUM comparators used in the appendix. The appendix power and ADL comparisons should therefore be interpreted as conditional power-delay comparisons rather than as comparisons among equally calibrated procedures."
  )
}

make_alt_metric_comment <- function(power_dt, scenario_value, metric) {
  sub <- power_dt[scenario == scenario_value & delta > 0]
  if (!nrow(sub)) return("")
  lead_det <- top_metric_row(sub, metric, c("detector"))
  lag_det <- bottom_metric_row(sub, metric, c("detector"))
  if (is.null(lead_det) || is.null(lag_det)) return("")
  det_phrase <- if (identical(metric, "adl")) {
    dominant_winner_phrase(winner_by_t(sub, metric, c("detector")), "detector", detector_pretty, verb = "has the smallest delay")
  } else {
    dominant_winner_phrase(winner_by_t(sub, metric, c("detector")), "detector", detector_pretty, verb = "leads")
  }
  progress <- progress_sentence(
    metric_path_t(sub[detector == lead_det$detector], metric),
    metric,
    detector_pretty[[lead_det$detector]],
    scenario_value
  )
  if (identical(metric, "adl")) {
    paste(
      "The retained alternative-detector ADL tables and subfigures favor the more aggressive RSMS benchmark shapes.",
      sprintf("Across detector families, %s.", det_phrase),
      progress,
      sprintf(
        "%s is the slowest retained detector family on average. Consistent with the main KS/CvM evidence, the SSMS rows mainly serve as conservative controls rather than as competitive supplementary detectors.",
        detector_pretty[[lag_det$detector]]
      ),
      "The current public replication bundle retains the Page-CUSUM, weighted-CUSUM, and multiscale-MOSUM comparators used in the appendix."
    )
  } else {
    paste(
      sprintf(
        "The retained alternative-detector %s tables and subfigures are led by %s.",
        if (identical(metric, "power")) "raw-power" else "size-adjusted-power",
        detector_pretty[[lead_det$detector]]
      ),
      sprintf("Across detector families, %s.", det_phrase),
      progress,
      "Across these retained detector shapes, weighted CUSUM and multiscale MOSUM are the empirically most informative supplementary comparators, whereas Page-CUSUM remains useful but typically secondary to those two.",
      "The corresponding SSMS rows remain materially more conservative than their RSMS counterparts, so they are best read as conservative controls rather than as competitive supplementary detectors.",
      sprintf(
        "%s is the weakest retained detector family on this criterion.",
        detector_pretty[[lag_det$detector]]
      )
    )
  }
}

main_ks_file <- function(metric, scenario_value, dgp_value) {
  sprintf("sim_%s_ks_%s_%s.png", scenario_value, metric, dgp_slug(dgp_value))
}

main_cvm_file <- function(metric, scenario_value, weight_value, dgp_value) {
  sprintf("sim_%s_cvm_%s_%s_%s.png", scenario_value, tolower(weight_value), metric, dgp_slug(dgp_value))
}

alt_file <- function(metric, scenario_value, tuning_slug, dgp_value) {
  sprintf("sim_alt_%s_%s_%s_%s.png", scenario_value, tuning_slug, metric, dgp_slug(dgp_value))
}

tuning_slug <- function(x) {
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  y <- gsub("^_+|_+$", "", y)
  y
}

null_files <- sort(list.files(null_dir, pattern = "\\.csv$", full.names = TRUE))
power_files <- sort(list.files(power_root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE))

main_null <- summarize_null_files(null_files, filter_main_dt)
alt_null <- summarize_null_files(null_files, filter_alt_dt)

if (nrow(main_null) == 0L) stop("No main-method null summaries were produced.")
if (nrow(alt_null) == 0L) stop("No alternative-detector null summaries were produced.")

main_meta <- build_method_meta(main_null, "Main")
alt_meta <- build_method_meta(alt_null, "Benchmark")

main_null <- merge(main_null, main_meta[, .(method_id, paper_label, panel, method_rank)], by = "method_id", all.x = TRUE, sort = FALSE)
alt_null <- merge(alt_null, alt_meta[, .(method_id, paper_label, panel, method_rank, tuning_label, tuning_label_plot, detector_rank)], by = "method_id", all.x = TRUE, sort = FALSE)

main_power <- summarize_power_files(power_files, filter_main_dt, main_null)
alt_power <- summarize_power_files(power_files, filter_alt_dt, alt_null)

if (nrow(main_power) == 0L) stop("No main-method power summaries were produced.")
if (nrow(alt_power) == 0L) stop("No alternative-detector power summaries were produced.")

main_power <- merge(main_power, main_meta[, .(method_id, paper_label, panel, method_rank)], by = "method_id", all.x = TRUE, sort = FALSE)
alt_power <- merge(alt_power, alt_meta[, .(method_id, paper_label, panel, method_rank, tuning_label, tuning_label_plot, detector_rank)], by = "method_id", all.x = TRUE, sort = FALSE)

fwrite(main_null, file.path(paper_summary_dir, "main_null_summary_alpha_005.csv"))
fwrite(main_power, file.path(paper_summary_dir, "main_power_curve_summary_alpha_005.csv"))
fwrite(alt_null, file.path(paper_summary_dir, "alt_null_summary_alpha_005.csv"))
fwrite(alt_power, file.path(paper_summary_dir, "alt_power_curve_summary_alpha_005.csv"))

best_selfnorm_power <- aggregate_metric_dt(
  main_power[scenario %in% main_scenarios & standardizer %in% c("RSMS", "SSMS") & delta > 0 & type == "KS"],
  metric = "sap",
  by_cols = c("paper_label")
)[order(-value)]

best_selfnorm_adl <- aggregate_metric_dt(
  main_power[scenario %in% main_scenarios & standardizer %in% c("RSMS", "SSMS") & delta > 0 & type == "KS"],
  metric = "adl",
  by_cols = c("paper_label")
)[order(value)]

best_selfnorm_power_label <- if (nrow(best_selfnorm_power) >= 1L) best_selfnorm_power$paper_label[1L] else "RSMS"
best_selfnorm_adl_label <- if (nrow(best_selfnorm_adl) >= 1L) best_selfnorm_adl$paper_label[1L] else "RSMS"

main_null_ks_table <- null_detail_table_lines(
  main_null,
  type_value = "KS",
  caption = "Null rejection rates (in \\%) at the 5\\% level for the KS monitor family with $m=500$. Columns keep the DGP and monitoring horizon $T\\in\\{1,2,5,10\\}$ explicit, so the remaining size distortions can be read design by design.",
  label = "tab:sim-null-ks",
  m_value = 500L
)

main_null_cvm_table <- null_detail_table_lines_cvm(
  main_null,
  caption = paste(
    "Null rejection rates (in \\%) at the 5\\% level for the weighted-CvM monitor family with $m=500$.",
    "Columns keep the DGP and monitoring horizon $T\\in\\{1,2,5,10\\}$ explicit across the four integrable weight families.",
    cvm_grouped_table_note
  ),
  label = "tab:sim-null-cvm",
  m_value = 500L
)

alt_null_table <- null_detail_table_core_lines(
  sub = alt_null,
  caption = "Null rejection rates (in \\%) at the 5\\% level for the FPCA-compressed Page-CUSUM, weighted-CUSUM, and multiscale MOSUM detector family with $m=500$. Columns keep the DGP and monitoring horizon $T\\in\\{1,2,5,10\\}$ explicit.",
  label = "tab:app-alt-null-detail",
  m_value = 500L
)

metrics <- c("power", "sap", "adl")
metric_heading <- c(
  power = "Raw power",
  sap = "Size-adjusted power",
  adl = "Average detection lag"
)

for (scenario_value in all_scenarios) {
  for (metric in metrics) {
    for (dgp_value in dgp_levels) {
      save_ks_dgp_plot(
        main_power,
        scenario_value = scenario_value,
        dgp_value = dgp_value,
        metric = metric,
        filename = main_ks_file(metric, scenario_value, dgp_value)
      )
    }
  }
}

for (scenario_value in all_scenarios) {
  for (weight_value in weight_levels) {
    for (metric in metrics) {
      for (dgp_value in dgp_levels) {
        save_cvm_dgp_plot(
          main_power,
          scenario_value = scenario_value,
          weight_value = weight_value,
          dgp_value = dgp_value,
          metric = metric,
          filename = main_cvm_file(metric, scenario_value, weight_value, dgp_value)
        )
      }
    }
  }
}

alt_tuning_meta <- unique(alt_meta[order(detector_rank), .(tuning_label_plot, tuning_label)])
alt_tuning_order <- alt_tuning_meta$tuning_label_plot
alt_tuning_caption_map <- setNames(alt_tuning_meta$tuning_label, alt_tuning_meta$tuning_label_plot)

for (scenario_value in all_scenarios) {
  for (tuning_value in alt_tuning_order) {
    for (metric in metrics) {
      for (dgp_value in dgp_levels) {
        save_alt_dgp_plot(
          alt_power,
          scenario_value = scenario_value,
          tuning_value = tuning_value,
          dgp_value = dgp_value,
          metric = metric,
          filename = alt_file(metric, scenario_value, tuning_slug(tuning_value), dgp_value)
        )
      }
    }
  }
}

main_block <- c(
  "",
  "\\section{Simulation Studies\\label{sec:simulation}}",
  "",
  "This section reports design-level rejection probabilities and detection-delay summaries computed from the full Monte Carlo output rather than from endpoint summaries. Size-adjusted power is obtained by recalibrating each design with the corresponding empirical 5\\% null rejection rate whenever that comparison is more informative than the raw rejection rate.",
  "",
  "The non-null settings considered here are a level shift, a smooth change, and an abrupt localized change. Because a comparable sinusoidal-change design is unavailable for the present calibration, we restrict attention to the designs that can be evaluated on a common basis. The main text emphasizes the level-shift and smooth-change cases. To keep the body focused, the weight-by-weight weighted-CvM curves, the full detection-delay breakdowns for these two settings, the localized-change design, and the FPCA-compressed alternative detector family are collected in Appendix~\\ref{app:additional-geometries}.",
  "",
  "All power and ADL summaries below are computed with $m=500$, monitoring horizons $T\\in\\{1,2,5,10\\}$, break locations $s^{\\star}\\in\\{50,200\\}$, and 1000 Monte Carlo replications per design. The main-method family consists of HAC, SSMS, and RSMS KS monitors for both boundary exponents $\\gamma\\in\\{0,0.15\\}$ together with the four weighted-CvM rules under the uniform, early, mid, and late weights. The tables keep the DGP, the monitoring horizon, and the positive break magnitude $\\Delta$ explicit, averaging only over the two break locations. The plots are organized DGP by DGP, and within each subfigure the internal columns correspond to $T=1,2,5,10$.",
  "",
  "\\subsection{Null calibration}",
  "",
  "Tables~\\ref{tab:sim-null-ks} and \\ref{tab:sim-null-cvm} report null rejection rates design by design. The basic calibration ranking is clear: HAC is the most aggressive benchmark, SSMS is typically the tightest null-calibrated self-normalized rule, and RSMS is the more competitive self-normalized rule once power and detection delay are read together. That is the comparison that matters for the detailed break-level tables and curves below.",
  "",
  main_null_ks_table,
  "",
  main_null_ks_comment(main_null),
  "",
  main_null_cvm_table,
  "",
  main_null_cvm_comment(main_null),
  "",
  "\\subsection{Power and size-adjusted power}",
  "",
  "The detailed tables now keep the exact positive break grid within each DGP instead of collapsing the design to weak and strong endpoints. This makes the monotone signal-strength pattern explicit: within any fixed method and monitoring horizon, moving from left to right across the columns increases the break magnitude, and the rejection probability rises accordingly in both the raw and size-adjusted summaries. The same monotone pattern is visible in the grayscale subfigures, which are now generated DGP by DGP and then assembled with explicit subcaptions.",
  "",
  sprintf("Within the self-normalized class, the strongest average size-adjusted KS performance in the main designs considered here is delivered by %s. This is precisely where RSMS becomes attractive in practice: its size-adjusted rejection paths usually stay above the corresponding SSMS paths over a broad part of the break grid while avoiding the larger null distortions of HAC.", best_selfnorm_power_label),
  "",
  "The weighted-CvM tables reported below show the same message after re-tabulation. Once the break magnitude is held fixed within each DGP and horizon, the late and uniform weights are often the most responsive, but the RSMS versions remain the stronger self-normalized choices in the most relevant parts of the design. The appendix then records the full weight-by-weight weighted-CvM curves and the complete ADL decomposition.",
  ""
)

for (scenario_value in main_scenarios) {
  main_block <- c(
    main_block,
    sprintf("\\subsubsection{%s}", scenario_pretty[[scenario_value]]),
    "",
    sprintf("Under the %s setting, the raw and size-adjusted summaries now keep the DGP, the monitoring horizon $T$, and the positive break magnitude $\\Delta$ explicit. The raw and recalibrated rejection paths both rise as the break magnitude increases, while the detailed subfigures make clear where RSMS stays above SSMS without inheriting the larger null distortions of HAC.", scenario_pretty_adj[[scenario_value]]),
    "",
    "\\medskip\\noindent\\textbf{Raw power.}",
    "",
    design_sidewaystable_lines(
      main_design_table_dt(main_power, scenario_value, "KS", "power"),
      caption = metric_design_table_caption("power", "the KS monitors", scenario_value),
      label = sprintf("tab:%s-ks-power", scenario_value),
      value_type = "pct"
    ),
    "",
    three_dgp_figure_lines(
      paths = fig_rel(vapply(dgp_levels, function(dg) main_ks_file("power", scenario_value, dg), character(1L))),
      figure_caption = sprintf(
        "%s for the KS monitors under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
        metric_curve_caption("power"),
        scenario_pretty_adj[[scenario_value]],
        ks_style_caption
      ),
      label = sprintf("fig:%s-ks-power", scenario_value)
    ),
    "",
    make_ks_metric_comment(main_power, scenario_value, "power"),
    "",
    design_sidewaystable_lines_cvm(
      cvm_design_table_dt(main_power, scenario_value, "power"),
      caption = paste(
        metric_design_table_caption("power", "the weighted-CvM monitors across all four weight families", scenario_value),
        cvm_grouped_table_note
      ),
      label = sprintf("tab:%s-cvm-power", scenario_value),
      value_type = "pct"
    ),
    "",
    make_cvm_metric_comment(main_power, scenario_value, "power"),
    "",
    sprintf("Appendix~\\ref{app:main-scenario-power} reports the weight-by-weight weighted-CvM raw-power curves for the %s setting.", scenario_pretty_adj[[scenario_value]]),
    "",
    "\\medskip\\noindent\\textbf{Size-adjusted power.}",
    "",
    design_sidewaystable_lines(
      main_design_table_dt(main_power, scenario_value, "KS", "sap"),
      caption = metric_design_table_caption("sap", "the KS monitors", scenario_value),
      label = sprintf("tab:%s-ks-sap", scenario_value),
      value_type = "pct"
    ),
    "",
    three_dgp_figure_lines(
      paths = fig_rel(vapply(dgp_levels, function(dg) main_ks_file("sap", scenario_value, dg), character(1L))),
      figure_caption = sprintf(
        "%s for the KS monitors under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
        metric_curve_caption("sap"),
        scenario_pretty_adj[[scenario_value]],
        ks_style_caption
      ),
      label = sprintf("fig:%s-ks-sap", scenario_value)
    ),
    "",
    make_ks_metric_comment(main_power, scenario_value, "sap"),
    "",
    design_sidewaystable_lines_cvm(
      cvm_design_table_dt(main_power, scenario_value, "sap"),
      caption = paste(
        metric_design_table_caption("sap", "the weighted-CvM monitors across all four weight families", scenario_value),
        cvm_grouped_table_note
      ),
      label = sprintf("tab:%s-cvm-sap", scenario_value),
      value_type = "pct"
    ),
    "",
    make_cvm_metric_comment(main_power, scenario_value, "sap"),
    "",
    sprintf("Appendix~\\ref{app:main-scenario-power} records the corresponding weight-by-weight weighted-CvM size-adjusted-power curves for the %s setting.", scenario_pretty_adj[[scenario_value]]),
    "",
    "\\FloatBarrier",
    ""
  )
}

main_block <- c(
  main_block,
  "\\subsection{Detection delay summary}",
  "",
  sprintf("To keep the body concise, the full ADL tables and curves for the level-shift and smooth-change settings are reported in Appendix~\\ref{app:main-scenario-adl}. The central implication is stable across those designs: within the self-normalized KS class, %s delivers the shortest delay over the broadest part of the break grid, whereas within the weighted-CvM family the earlier weights respond faster but do not match the raw- and size-adjusted-power gains delivered by the late and uniform weights.", best_selfnorm_adl_label),
  "",
  "This presentation leaves the core calibration and power comparison in the main text while preserving the complete timing decomposition for readers who wish to inspect the break-level delay patterns in detail.",
  ""
)

appendix_block <- c(
  "",
  "\\section{Additional simulation evidence}\\label{app:additional-geometries}",
  "",
  "This appendix reports three additional layers of simulation evidence. First, it collects the supplementary weight-by-weight weighted-CvM power curves and the full ADL breakdowns for the level-shift and smooth-change settings emphasized in the main text. Second, it reports the abrupt localized-change setting for the main HAC, SSMS, and RSMS monitors. Third, it tabulates and plots a retained set of FPCA-compressed alternative detector shapes built from weighted CUSUM, Page-CUSUM, and multiscale MOSUM statistics computed from the same retained score vectors.",
  "",
  "\\subsection{Supplementary power evidence for the level-shift and smooth-change settings}\\label{app:main-scenario-power}",
  "",
  "The main text reports the consolidated KS figures and the combined weighted-CvM tables for the two principal scenarios. This subsection records the weight-by-weight weighted-CvM raw-power and size-adjusted-power curves so that the full graphical evidence remains available without lengthening the body of the paper.",
  ""
)

for (scenario_value in main_scenarios) {
  appendix_block <- c(
    appendix_block,
    sprintf("\\subsubsection{%s}", scenario_pretty[[scenario_value]]),
    "",
    "\\medskip\\noindent\\textbf{Raw power.}",
    ""
  )
  for (weight_value in weight_levels) {
    appendix_block <- c(
      appendix_block,
      three_dgp_figure_lines(
        paths = fig_rel(vapply(dgp_levels, function(dg) main_cvm_file("power", scenario_value, weight_value, dg), character(1L))),
        figure_caption = sprintf(
          "%s for the weighted-CvM monitors using %s under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
          metric_curve_caption("power"),
          weight_caption_pretty[[weight_value]],
          scenario_pretty_adj[[scenario_value]],
          std_style_caption
        ),
        label = sprintf("fig:%s-cvm-%s-power", scenario_value, tolower(weight_value))
      ),
      ""
    )
  }
  appendix_block <- c(
    appendix_block,
    make_cvm_metric_comment(main_power, scenario_value, "power"),
    "",
    "\\medskip\\noindent\\textbf{Size-adjusted power.}",
    ""
  )
  for (weight_value in weight_levels) {
    appendix_block <- c(
      appendix_block,
      three_dgp_figure_lines(
        paths = fig_rel(vapply(dgp_levels, function(dg) main_cvm_file("sap", scenario_value, weight_value, dg), character(1L))),
        figure_caption = sprintf(
          "%s for the weighted-CvM monitors using %s under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
          metric_curve_caption("sap"),
          weight_caption_pretty[[weight_value]],
          scenario_pretty_adj[[scenario_value]],
          std_style_caption
        ),
        label = sprintf("fig:%s-cvm-%s-sap", scenario_value, tolower(weight_value))
      ),
      ""
    )
  }
  appendix_block <- c(
    appendix_block,
    make_cvm_metric_comment(main_power, scenario_value, "sap"),
    "",
    "\\FloatBarrier",
    ""
  )
}

appendix_block <- c(
  appendix_block,
  "\\subsection{Average detection lag for the level-shift and smooth-change settings}\\label{app:main-scenario-adl}",
  "",
  "This subsection collects the full ADL breakdown for the two scenarios emphasized in the main text. The tables keep the break grid explicit, and the accompanying grayscale subfigures preserve the DGP-by-DGP layout used throughout the simulation section.",
  "",
  ""
)

for (scenario_value in main_scenarios) {
  appendix_block <- c(
    appendix_block,
    sprintf("\\subsubsection{%s}", scenario_pretty[[scenario_value]]),
    "",
    design_sidewaystable_lines(
      main_design_table_dt(main_power, scenario_value, "KS", "adl"),
      caption = metric_design_table_caption("adl", "the KS monitors", scenario_value),
      label = sprintf("tab:%s-ks-adl", scenario_value),
      value_type = "num"
    ),
    "",
    three_dgp_figure_lines(
      paths = fig_rel(vapply(dgp_levels, function(dg) main_ks_file("adl", scenario_value, dg), character(1L))),
      figure_caption = sprintf(
        "%s for the KS monitors under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
        metric_curve_caption("adl"),
        scenario_pretty_adj[[scenario_value]],
        ks_style_caption
      ),
      label = sprintf("fig:%s-ks-adl", scenario_value)
    ),
    "",
    make_ks_metric_comment(main_power, scenario_value, "adl"),
    "",
    design_sidewaystable_lines_cvm(
      cvm_design_table_dt(main_power, scenario_value, "adl"),
      caption = paste(
        metric_design_table_caption("adl", "the weighted-CvM monitors across all four weight families", scenario_value),
        cvm_grouped_table_note
      ),
      label = sprintf("tab:%s-cvm-adl", scenario_value),
      value_type = "num"
    ),
    ""
  )
  for (weight_value in weight_levels) {
    appendix_block <- c(
      appendix_block,
      three_dgp_figure_lines(
        paths = fig_rel(vapply(dgp_levels, function(dg) main_cvm_file("adl", scenario_value, weight_value, dg), character(1L))),
        figure_caption = sprintf(
          "%s for the weighted-CvM monitors using %s under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
          metric_curve_caption("adl"),
          weight_caption_pretty[[weight_value]],
          scenario_pretty_adj[[scenario_value]],
          std_style_caption
        ),
        label = sprintf("fig:%s-cvm-%s-adl", scenario_value, tolower(weight_value))
      ),
      ""
    )
  }
  appendix_block <- c(
    appendix_block,
    make_cvm_metric_comment(main_power, scenario_value, "adl"),
    "",
    "\\FloatBarrier",
    ""
  )
}

appendix_block <- c(
  appendix_block,
  "\\subsection{Localized break setting for the main monitors}",
  "",
  "The localized-change appendix keeps the same break-level organization as the main text. The tables retain the exact positive break magnitudes within each DGP and horizon, while the grayscale subfigures are again generated DGP by DGP so that each subcaption is unambiguous.",
  ""
)

for (metric in metrics) {
  appendix_block <- c(
    appendix_block,
    sprintf("\\subsubsection{%s}", metric_heading[[metric]]),
    "",
    design_sidewaystable_lines(
      main_design_table_dt(main_power, "abrupt_local_change", "KS", metric),
      caption = metric_design_table_caption(metric, "the KS monitors", "abrupt_local_change"),
      label = sprintf("tab:abrupt-ks-%s", metric),
      value_type = if (metric %in% c("power", "sap")) "pct" else "num"
    ),
    "",
    three_dgp_figure_lines(
      paths = fig_rel(vapply(dgp_levels, function(dg) main_ks_file(metric, "abrupt_local_change", dg), character(1L))),
      figure_caption = sprintf(
        "%s for the KS monitors under the localized-change setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
        metric_curve_caption(metric),
        ks_style_caption
      ),
      label = sprintf("fig:abrupt-ks-%s", metric)
    ),
    "",
    make_ks_metric_comment(main_power, "abrupt_local_change", metric),
    "",
    design_sidewaystable_lines_cvm(
      cvm_design_table_dt(main_power, "abrupt_local_change", metric),
      caption = paste(
        metric_design_table_caption(metric, "the weighted-CvM monitors across all four weight families", "abrupt_local_change"),
        cvm_grouped_table_note
      ),
      label = sprintf("tab:abrupt-cvm-%s", metric),
      value_type = if (metric %in% c("power", "sap")) "pct" else "num"
    ),
    ""
  )
  for (weight_value in weight_levels) {
    appendix_block <- c(
      appendix_block,
      three_dgp_figure_lines(
        paths = fig_rel(vapply(dgp_levels, function(dg) main_cvm_file(metric, "abrupt_local_change", weight_value, dg), character(1L))),
        figure_caption = sprintf(
          "%s for the weighted-CvM monitors using %s under the localized-change setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
          metric_curve_caption(metric),
          weight_caption_pretty[[weight_value]],
          std_style_caption
        ),
        label = sprintf("fig:abrupt-cvm-%s-%s", tolower(weight_value), metric)
      ),
      ""
    )
  }
  appendix_block <- c(
    appendix_block,
    make_cvm_metric_comment(main_power, "abrupt_local_change", metric),
    "",
    "\\FloatBarrier",
    ""
  )
}

appendix_block <- c(
  appendix_block,
  "\\subsection{Alternative FPCA-compressed detector family}",
  "",
  "Table~\\ref{tab:app-alt-null-detail} reports the null calibration of the retained alternative detector family. The scenario-specific tables and subfigures that follow then keep the raw DGP-specific break grids visible, so the weighted-CUSUM, Page-CUSUM, and multiscale MOSUM comparisons are no longer compressed to weak and strong endpoints.",
  "",
  alt_null_table,
  "",
  alt_null_comment(alt_null),
  ""
)

for (scenario_value in all_scenarios) {
  appendix_block <- c(
    appendix_block,
    sprintf("\\subsubsection{%s}", scenario_pretty[[scenario_value]]),
    "",
    sprintf("Under the %s setting, the appendix tables retain the DGP, the horizon $T$, and the positive break grid for the FPCA-compressed weighted-CUSUM, Page-CUSUM, and multiscale MOSUM detectors. The accompanying grayscale figures are again generated DGP by DGP, with the four internal columns corresponding to $T=1,2,5,10$.", scenario_pretty_adj[[scenario_value]]),
    ""
  )
  for (metric in metrics) {
    appendix_block <- c(
      appendix_block,
      sprintf("\\medskip\\noindent\\textbf{%s.}", metric_heading[[metric]]),
      "",
      design_sidewaystable_lines(
        alt_design_table_dt(alt_power, scenario_value, metric),
        caption = metric_design_table_caption(metric, "the FPCA-compressed alternative detector family", scenario_value),
        label = sprintf("tab:alt-%s-%s", scenario_value, metric),
        value_type = if (metric %in% c("power", "sap")) "pct" else "num"
      ),
      ""
    )
    for (tuning_value in alt_tuning_order) {
      appendix_block <- c(
        appendix_block,
        three_dgp_figure_lines(
          paths = fig_rel(vapply(dgp_levels, function(dg) alt_file(metric, scenario_value, tuning_slug(tuning_value), dg), character(1L))),
          figure_caption = sprintf(
            "%s for %s under the %s setting. Each subfigure fixes one DGP, and the four internal columns correspond to $T=1,2,5,10$. %s",
            metric_curve_caption(metric),
            alt_tuning_caption_map[[tuning_value]],
            scenario_pretty_adj[[scenario_value]],
            std_style_caption
          ),
          label = sprintf("fig:alt-%s-%s-%s", scenario_value, tuning_slug(tuning_value), metric)
        ),
        ""
      )
    }
    appendix_block <- c(
      appendix_block,
      make_alt_metric_comment(alt_power, scenario_value, metric),
      ""
    )
  }
  appendix_block <- c(appendix_block, "\\FloatBarrier", "")
}

tex_lines <- readLines(tex_path, warn = FALSE)
tex_lines <- replace_between_markers(tex_lines, marker_main_start, marker_main_end, main_block)
tex_lines <- replace_between_markers(tex_lines, marker_app_start, marker_app_end, appendix_block)
writeLines(tex_lines, tex_path, useBytes = TRUE)

message("Refreshed simulation LaTeX injected into main manuscript.")
message("Summary CSVs written to: ", normalizePath(paper_summary_dir, winslash = "/", mustWork = FALSE))
message("Figures written to: ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE))
