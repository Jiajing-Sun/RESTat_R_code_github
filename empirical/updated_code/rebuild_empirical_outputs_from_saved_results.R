get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

na_blank <- function(x, blank = "") {
  x[is.na(x) | x == "NA"] <- blank
  x
}

fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  out <- ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
  unname(out)
}

fmt_yes_no <- function(x) {
  x <- toupper(as.character(x))
  ifelse(x == "TRUE", "Yes", "No")
}

fmt_signal_date <- function(x) {
  x <- na_blank(as.character(x))
  out <- rep("No signal", length(x))
  keep <- nzchar(x)
  out[keep] <- format(as.Date(x[keep]), "%d %b %Y")
  out
}

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

write_latex_table <- function(df, path, caption, label, align, size_cmd = "\\scriptsize", note = NULL) {
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    size_cmd,
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

build_rank_table <- function(df, method_col) {
  keep <- toupper(as.character(df$reject)) == "TRUE"
  out <- df[keep, c("T_label", method_col, "first_rejection_date", "trading_gap_from_2020_03_16", "stat_over_cv")]
  names(out)[2] <- "method_label"
  out$T_label <- as.integer(out$T_label)
  out$first_rejection_date <- as.Date(out$first_rejection_date)
  out$trading_gap_from_2020_03_16 <- as.integer(out$trading_gap_from_2020_03_16)
  out$stat_over_cv <- as.numeric(out$stat_over_cv)
  out <- out[order(out$T_label, out$first_rejection_date, -out$stat_over_cv), ]
  out$first_rejection_date <- format(out$first_rejection_date, "%Y-%m-%d")
  out
}

script_dir <- get_script_dir()
empirical_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(empirical_dir, "data")
out_dir <- file.path(empirical_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main_file <- file.path(data_dir, "empirical_main_results.csv")
alt_file <- file.path(data_dir, "empirical_alt_results.csv")

if (!file.exists(main_file) || !file.exists(alt_file)) {
  stop("Saved empirical result files are missing from ", data_dir)
}

main <- read.csv(main_file, stringsAsFactors = FALSE, check.names = FALSE)
alt <- read.csv(alt_file, stringsAsFactors = FALSE, check.names = FALSE)

main$T_label_num <- as.integer(main$T_label)
alt$T_label_num <- as.integer(alt$T_label)

std_order <- c("HAC", "RSMS", "SSMS")
weight_order <- c("U", "Early", "Mid", "Late")
detector_order <- c("WeightedCUSUM", "MultiscaleMOSUM", "PageCUSUM", "FullCUSUM", "MOSUM")

ks <- main[main$type == "KS", ]
ks <- ks[order(ks$T_label_num, match(ks$standardizer, std_order), as.numeric(ks$gamma)), ]
ks_out <- data.frame(
  T = ks$T_label,
  Method = ks$standardizer,
  Gamma = paste0("gamma=", ifelse(ks$gamma == "0", "0.00", ks$gamma)),
  Statistic = fmt_num(ks$statistic),
  `Critical value` = fmt_num(ks$critical_value),
  `Stat./CV` = fmt_num(ks$stat_over_cv),
  `Reject?` = fmt_yes_no(ks$reject),
  `First signal` = fmt_signal_date(ks$first_rejection_date),
  check.names = FALSE
)
write.csv(ks_out, file.path(out_dir, "empirical_main_ks_table.csv"), row.names = FALSE, na = "")
write_latex_table(
  ks_out,
  file.path(out_dir, "empirical_main_ks_table.tex"),
  "Empirical KS monitoring results rebuilt from the saved current-paper output files.",
  "tab:bundle-empirical-main-ks",
  "cccccccl"
)

cvm <- main[main$type == "CvM", ]
cvm <- cvm[order(cvm$T_label_num, match(cvm$weight_name, weight_order), match(cvm$standardizer, std_order)), ]
cvm_out <- data.frame(
  T = cvm$T_label,
  Weight = cvm$weight_name,
  Method = cvm$standardizer,
  Statistic = fmt_num(cvm$statistic),
  `Critical value` = fmt_num(cvm$critical_value),
  `Stat./CV` = fmt_num(cvm$stat_over_cv),
  `Reject?` = fmt_yes_no(cvm$reject),
  `First signal` = fmt_signal_date(cvm$first_rejection_date),
  check.names = FALSE
)
write.csv(cvm_out, file.path(out_dir, "empirical_main_cvm_table.csv"), row.names = FALSE, na = "")
write_latex_table(
  cvm_out,
  file.path(out_dir, "empirical_main_cvm_table.tex"),
  "Empirical weighted-CvM monitoring results rebuilt from the saved current-paper output files.",
  "tab:bundle-empirical-main-cvm",
  "cccccccl"
)

detector_label <- function(detector, gamma, bandwidth_h, omega_name, hset_name) {
  gamma <- na_blank(as.character(gamma))
  bandwidth_h <- na_blank(as.character(bandwidth_h))
  omega_name <- na_blank(as.character(omega_name))
  hset_name <- na_blank(as.character(hset_name))

  if (detector == "WeightedCUSUM") {
    return("Weighted CUSUM")
  }
  if (detector == "MultiscaleMOSUM") {
    return("Multiscale MOSUM")
  }
  if (detector == "PageCUSUM") {
    return("Page-CUSUM")
  }
  if (detector == "FullCUSUM") {
    return("Full-CUSUM")
  }
  if (detector == "MOSUM") {
    return(paste0("MOSUM (h=", bandwidth_h, ")"))
  }
  detector
}

gamma_label <- function(gamma) {
  gamma <- na_blank(as.character(gamma))
  ifelse(
    nzchar(gamma),
    paste0("gamma=", ifelse(gamma == "0", "0.00", gamma)),
    ""
  )
}

alt$DetectorLabel <- mapply(
  detector_label,
  alt$detector,
  alt$gamma,
  alt$bandwidth_h,
  alt$omega_name,
  alt$hset_name,
  USE.NAMES = FALSE
)

alt_full <- alt[order(alt$T_label_num, match(alt$detector, detector_order), match(alt$standardizer, std_order), as.numeric(ifelse(alt$gamma == "NA", NA, alt$gamma))), ]
alt_full_out <- data.frame(
  T = alt_full$T_label,
  Detector = alt_full$DetectorLabel,
  Method = alt_full$standardizer,
  Gamma = gamma_label(alt_full$gamma),
  h = na_blank(alt_full$bandwidth_h),
  Omega = na_blank(alt_full$omega_name),
  Hset = na_blank(alt_full$hset_name),
  Statistic = fmt_num(alt_full$statistic),
  `Critical value` = fmt_num(alt_full$critical_value),
  `Stat./CV` = fmt_num(alt_full$stat_over_cv),
  `Reject?` = fmt_yes_no(alt_full$reject),
  `First signal` = fmt_signal_date(alt_full$first_rejection_date),
  check.names = FALSE
)
write.csv(alt_full_out, file.path(out_dir, "empirical_alt_full_table.csv"), row.names = FALSE, na = "")
write_latex_table(
  alt_full_out,
  file.path(out_dir, "empirical_alt_full_table.tex"),
  "Full empirical alternative-detector results rebuilt from the saved current-paper output files.",
  "tab:bundle-empirical-alt-full",
  "cccccccccccl",
  note = "The retained supplement discussion emphasizes weighted CUSUM and multiscale MOSUM, with Page-CUSUM kept as a secondary comparator."
)

keep_detectors <- c("WeightedCUSUM", "MultiscaleMOSUM", "PageCUSUM")
alt_retained <- alt[alt$detector %in% keep_detectors, ]
alt_retained <- alt_retained[order(alt_retained$T_label_num, match(alt_retained$detector, keep_detectors), match(alt_retained$standardizer, std_order), as.numeric(ifelse(alt_retained$gamma == "NA", NA, alt_retained$gamma))), ]
alt_retained_out <- data.frame(
  T = alt_retained$T_label,
  Detector = alt_retained$DetectorLabel,
  Method = alt_retained$standardizer,
  Gamma = gamma_label(alt_retained$gamma),
  h = na_blank(alt_retained$bandwidth_h),
  Omega = na_blank(alt_retained$omega_name),
  Hset = na_blank(alt_retained$hset_name),
  Statistic = fmt_num(alt_retained$statistic),
  `Critical value` = fmt_num(alt_retained$critical_value),
  `Stat./CV` = fmt_num(alt_retained$stat_over_cv),
  `Reject?` = fmt_yes_no(alt_retained$reject),
  `First signal` = fmt_signal_date(alt_retained$first_rejection_date),
  check.names = FALSE
)
write.csv(alt_retained_out, file.path(out_dir, "empirical_alt_retained_table.csv"), row.names = FALSE, na = "")
write_latex_table(
  alt_retained_out,
  file.path(out_dir, "empirical_alt_retained_table.tex"),
  "Retained empirical alternative-detector results for the current supplement emphasis.",
  "tab:bundle-empirical-alt-retained",
  "cccccccccccl",
  note = "This retained table keeps weighted CUSUM and multiscale MOSUM as the detailed benchmarks and Page-CUSUM as the secondary comparison."
)

main_rank <- build_rank_table(main, "method_label")
alt_rank <- build_rank_table(alt[alt$detector %in% keep_detectors, ], "method_label")
write.csv(main_rank, file.path(out_dir, "empirical_main_detection_rank.csv"), row.names = FALSE, na = "")
write.csv(alt_rank, file.path(out_dir, "empirical_alt_retained_detection_rank.csv"), row.names = FALSE, na = "")

note_path <- file.path(out_dir, "empirical_narrative_notes.txt")
training_start <- unique(main$training_start)
training_end <- unique(main$training_end)
monitoring_start <- unique(main$monitoring_start)
monitoring_end <- unique(main$monitoring_end)
q_used <- unique(main$q_used)
m_train <- unique(main$m_train)
t_values <- sort(unique(main$T_label_num))

main_top <- head(main_rank, 12)
alt_top <- head(alt_rank, 12)

notes <- c(
  paste0("Empirical design: m = ", m_train[1], ", q_used = ", q_used[1], ", finite horizons = ", paste(t_values, collapse = ", "), "."),
  paste0("Training window: ", training_start[1], " to ", training_end[1], "."),
  paste0("Monitoring window: ", monitoring_start[1], " to ", monitoring_end[1], "."),
  "",
  "Earliest main-family detections by horizon:",
  capture.output(print(main_top, row.names = FALSE)),
  "",
  "Earliest retained alternative-detector detections by horizon:",
  capture.output(print(alt_top, row.names = FALSE))
)
writeLines(notes, con = note_path, useBytes = TRUE)

message("Saved rebuilt empirical outputs to: ", out_dir)
