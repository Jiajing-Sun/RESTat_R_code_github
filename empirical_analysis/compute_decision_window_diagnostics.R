decode_rscript_path <- function(x) {
  gsub("~\\+~", " ", x)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(decode_rscript_path(sub("^--file=", "", file_arg[1])))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

write_latex_table_simple <- function(df, path, caption, label, align, note = NULL) {
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  if (!is.null(note)) {
    lines <- c(
      lines,
      "\\begin{flushleft}",
      "\\footnotesize",
      note,
      "\\end{flushleft}"
    )
  }
  lines <- c(lines, "\\end{table}")
  writeLines(lines, con = path, useBytes = TRUE)
}

fmt_num <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits)
}

script_dir <- get_script_dir()
empirical_dir <- normalizePath(script_dir)
data_file <- Sys.getenv("SPX_DATA_FILE", unset = file.path(empirical_dir, "data", "SPX.csv"))
out_dir <- file.path(empirical_dir, "outputs")

if (!file.exists(data_file)) {
  stop("Could not find SPX.csv at: ", data_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spx <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
blank_names <- names(spx) == ""
if (any(blank_names)) {
  spx <- spx[, !blank_names, drop = FALSE]
}
if (!("DateTime" %in% names(spx)) || !("Close" %in% names(spx))) {
  stop("SPX.csv must contain DateTime and Close columns.")
}

spx$DateTime <- as.POSIXct(spx$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
spx$Date <- as.Date(spx$DateTime)
spx$Close <- as.numeric(spx$Close)
spx <- spx[!is.na(spx$Date) & !is.na(spx$Close), ]
spx <- spx[spx$Date >= as.Date("2019-10-21") & spx$Date <= as.Date("2020-12-31"), ]

day_split <- split(spx$Close, spx$Date)
daily_metrics <- do.call(
  rbind,
  lapply(names(day_split), function(day_name) {
    close_vec <- as.numeric(day_split[[day_name]])
    if (length(close_vec) < 2) {
      return(NULL)
    }
    ret_vec <- diff(log(close_vec))
    data.frame(
      Date = as.Date(day_name),
      avg_abs_1min = mean(abs(ret_vec), na.rm = TRUE),
      daily_realized_variation = sum(ret_vec^2, na.rm = TRUE),
      daily_log_return = sum(ret_vec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
daily_metrics <- daily_metrics[order(daily_metrics$Date), , drop = FALSE]

summarize_window <- function(start_date, end_date) {
  keep <- daily_metrics$Date >= as.Date(start_date) & daily_metrics$Date <= as.Date(end_date)
  window_df <- daily_metrics[keep, , drop = FALSE]
  c(
    trading_days = nrow(window_df),
    avg_abs_1min = mean(window_df$avg_abs_1min),
    avg_daily_rv = mean(window_df$daily_realized_variation),
    cum_log_return = sum(window_df$daily_log_return)
  )
}

window_specs <- list(
  list(
    label = "Training sample",
    start = "2019-10-21",
    end = "2019-12-31"
  ),
  list(
    label = "Early 2020 before first main alarm",
    start = "2020-01-02",
    end = "2020-02-26"
  ),
  list(
    label = "RSMS(T=2) to HAC(T=2)",
    start = "2020-02-27",
    end = "2020-04-01"
  ),
  list(
    label = "RSMS(T=5) to HAC(T=5)",
    start = "2020-03-13",
    end = "2020-04-02"
  )
)

window_summary <- do.call(
  rbind,
  lapply(window_specs, function(spec) {
    stats <- summarize_window(spec$start, spec$end)
    data.frame(
      window = spec$label,
      start_date = spec$start,
      end_date = spec$end,
      trading_days = unname(stats["trading_days"]),
      avg_abs_1min = unname(stats["avg_abs_1min"]),
      avg_daily_rv = unname(stats["avg_daily_rv"]),
      cum_log_return = unname(stats["cum_log_return"]),
      stringsAsFactors = FALSE
    )
  })
)

training_abs <- window_summary$avg_abs_1min[window_summary$window == "Training sample"][1]
training_rv <- window_summary$avg_daily_rv[window_summary$window == "Training sample"][1]

window_summary$avg_abs_1min_bp <- window_summary$avg_abs_1min * 10000
window_summary$avg_daily_rv_x1e4 <- window_summary$avg_daily_rv * 10000
window_summary$rel_abs_to_training <- window_summary$avg_abs_1min / training_abs
window_summary$rel_rv_to_training <- window_summary$avg_daily_rv / training_rv

write.csv(
  window_summary,
  file.path(out_dir, "empirical_decision_window_metrics.csv"),
  row.names = FALSE,
  na = ""
)

latex_df <- data.frame(
  Window = c(
    "Training sample",
    "Early 2020 pre-alarm period",
    "RSMS($T=2$) to HAC($T=2$)",
    "RSMS($T=5$) to HAC($T=5$)"
  ),
  `Trading days` = as.integer(window_summary$trading_days),
  `Avg. $|r|$ (bp)` = fmt_num(window_summary$avg_abs_1min_bp, 2),
  `Avg. daily RV ($\\times 10^{-4}$)` = fmt_num(window_summary$avg_daily_rv_x1e4, 2),
  `Rel. $|r|$` = paste0(fmt_num(window_summary$rel_abs_to_training, 1), "$\\times$"),
  `Rel. RV` = paste0(fmt_num(window_summary$rel_rv_to_training, 1), "$\\times$"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  latex_df,
  file.path(out_dir, "empirical_decision_window_table.csv"),
  row.names = FALSE,
  na = ""
)

write_latex_table_simple(
  latex_df,
  file.path(out_dir, "empirical_decision_window_table.tex"),
  "Raw intraday stress diagnostics for the main KS decision windows.",
  "tab:bundle-empirical-decision-window",
  "lccccc",
  note = paste(
    "The diagnostics are computed directly from the raw 1-minute S\\&P 500 series.",
    "Average absolute returns are reported in basis points and daily realized variation is scaled by $10^{-4}$.",
    "The four windows are 2019-10-21 to 2019-12-31, 2020-01-02 to 2020-02-26,",
    "2020-02-27 to 2020-04-01, and 2020-03-13 to 2020-04-02."
  )
)

message("Saved decision-window diagnostics to: ", out_dir)
