decode_rscript_path <- function(x) {
  gsub("~\\+~", " ", x)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(decode_rscript_path(sub("^--file=", "", file_arg[1])))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

fmt_signal_date <- function(x) {
  x <- as.character(x)
  out <- rep("No signal", length(x))
  keep <- !is.na(x) & nzchar(x)
  if (any(keep)) {
    d <- as.Date(x[keep])
    out[keep] <- paste0(format(d, "%b."), "~", as.integer(format(d, "%d")), ", ", format(d, "%Y"))
  }
  out
}

write_latex_table <- function(df, path) {
  align <- paste0("l", paste(rep("c", ncol(df) - 1L), collapse = ""))
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.6pt}",
    "\\renewcommand{\\arraystretch}{0.95}",
    "\\caption{Alternative training-window robustness for the main empirical KS alarms.}",
    "\\label{tab:empirical-training-window-robustness}",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{flushleft}",
    "\\footnotesize",
    "Note: Each row re-estimates the training-window FPCA basis using the last $m$ pre-2020 trading days ending on December~31, 2019, retains the number of components selected by the 80\\% FVE rule, and reports first 5\\% KS alarm dates. SSMS remains below the threshold in all rows. The check is a robustness exercise for the empirical ranking, not a replacement for the $m=50$, $K=4$ baseline used in the paper.",
    "\\end{flushleft}",
    "\\end{table}"
  )
  writeLines(lines, con = path, useBytes = TRUE)
}

build_fd_data_for_training_window <- function(spx, m, max_T = 5L,
                                              training_end = as.Date("2019-12-31")) {
  all_dates <- sort(unique(spx$Date))
  train_pool <- all_dates[all_dates <= training_end]
  if (length(train_pool) < m) {
    stop("Not enough pre-2020 trading days for m=", m)
  }
  training_dates <- tail(train_pool, m)
  monitor_pool <- all_dates[all_dates > training_end]
  needed_monitor <- m * max_T
  if (length(monitor_pool) < needed_monitor) {
    stop("Not enough monitoring days for m=", m, ", T=", max_T)
  }
  selected_dates <- c(training_dates, head(monitor_pool, needed_monitor))
  sub <- spx[spx$Date %in% selected_dates, c("Date", "DateTime", "Close")]
  sub <- sub[order(sub$Date, sub$DateTime), ]

  daily_list <- lapply(selected_dates, function(d) {
    day <- sub[sub$Date == d, c("Date", "DateTime", "Close")]
    day$Close <- c(NA_real_, diff(log(day$Close)))
    day <- day[!is.na(day$Close), ]
    rownames(day) <- NULL
    day
  })

  min_obs <- min(vapply(daily_list, nrow, integer(1L)))
  nbasis <- min(21L, floor(min_obs * 0.8))
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), norder = 4, nbasis = nbasis)
  coefs <- matrix(nrow = nbasis, ncol = length(daily_list))

  for (i in seq_along(daily_list)) {
    day <- daily_list[[i]]
    fd_obj <- fda::smooth.basis(
      argvals = seq(0, 1, length.out = nrow(day)),
      y = day$Close,
      fdParobj = basis
    )
    coefs[, i] <- fd_obj$fd$coefs
  }

  list(
    fd_data = fda::fd(
      coef = coefs,
      basisobj = basis,
      fdnames = list(time = "t", replicates = as.character(selected_dates), values = "SPX_Close")
    ),
    dates = selected_dates
  )
}

evaluate_training_window <- function(fd_data, all_dates, cv_main, m,
                                     T_grid = c(1L, 2L, 5L),
                                     gamma_grid = c(0, 0.15),
                                     alpha = 0.05) {
  main_catalog <- build_main_method_catalog(gamma_vec = gamma_grid)
  main_catalog <- subset(
    main_catalog,
    type == "KS" &
      standardizer %in% c("HAC", "SSMS", "RSMS") &
      gamma %in% gamma_grid
  )

  fp <- fpca_project_scores(fd_data, m = m, fixed_q = NA_integer_, fve_threshold = 0.80, n_grid = 301L)
  white <- whiten_scores(fp$scores_train, fp$scores_monitor)

  out <- list()
  idx <- 1L
  for (T_value in T_grid) {
    k_monitor <- m * T_value
    raw_train <- fp$scores_train
    raw_monitor <- fp$scores_monitor[seq_len(k_monitor), , drop = FALSE]
    rsms_train <- white$scores_train
    rsms_monitor <- white$scores_monitor[seq_len(k_monitor), , drop = FALSE]

    context_map <- list(
      HAC = prepare_standardizer_context(raw_train, raw_monitor, "HAC"),
      SSMS = prepare_standardizer_context(raw_train, raw_monitor, "SSMS"),
      RSMS = prepare_standardizer_context(rsms_train, rsms_monitor, "RSMS")
    )

    for (i in seq_len(nrow(main_catalog))) {
      row <- main_catalog[i, , drop = FALSE]
      eval_out <- evaluate_main_method(
        row,
        context_map,
        q_used = fp$q_used,
        T_value = T_value,
        alpha_levels = alpha,
        cv_main = cv_main
      )
      eval_out$m <- m
      eval_out$q_used <- fp$q_used
      eval_out$T_value <- T_value
      eval_out$T_label <- canonical_T_scalar(T_value)
      eval_out$first_rejection_date <- ifelse(
        !is.na(eval_out$first_rejection),
        as.character(all_dates[m + eval_out$first_rejection]),
        NA
      )
      out[[idx]] <- eval_out
      idx <- idx + 1L
    }
  }
  do.call(rbind, out)
}

build_summary_table <- function(results) {
  sub <- results[results$standardizer %in% c("RSMS", "HAC"), , drop = FALSE]
  rows <- list()
  idx <- 1L
  for (m_value in sort(unique(sub$m))) {
    for (gamma_value in sort(unique(sub$gamma))) {
      hit <- sub[sub$m == m_value & sub$gamma == gamma_value, , drop = FALSE]
      q_used <- unique(hit$q_used)
      row <- data.frame(
        `$m$` = m_value,
        `$K$` = q_used[1L],
        `$\\gamma$` = ifelse(abs(gamma_value) < 1e-12, "0", "0.15"),
        `RSMS ($T=1$)` = "No signal",
        `RSMS ($T=2$)` = "No signal",
        `RSMS ($T=5$)` = "No signal",
        `HAC ($T=1$)` = "No signal",
        `HAC ($T=2$)` = "No signal",
        `HAC ($T=5$)` = "No signal",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      for (std in c("RSMS", "HAC")) {
        for (T_label in c("1", "2", "5")) {
          one <- hit[hit$standardizer == std & hit$T_label == T_label, , drop = FALSE]
          if (nrow(one) == 1L) {
            row[[paste0(std, " ($T=", T_label, "$)")]] <- fmt_signal_date(one$first_rejection_date)
          }
        }
      }
      rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

script_dir <- get_script_dir()
empirical_dir <- normalizePath(script_dir)
release_root <- normalizePath(file.path(empirical_dir, ".."))
main_modules_dir <- file.path(release_root, "simulations", "main", "R")

for (f in c("dependencies.R", "utils.R", "critical_values_lookup.R", "method_catalog.R",
            "fpca_pipeline.R", "detectors_main.R")) {
  source(file.path(main_modules_dir, f), local = TRUE)
}
ensure_simulation_packages(install_if_missing = FALSE)

spx_path <- Sys.getenv("SPX_DATA_FILE", unset = file.path(empirical_dir, "data", "SPX.csv"))
if (!file.exists(spx_path)) {
  stop("Could not find SPX.csv at: ", spx_path)
}
spx <- read.csv(spx_path, stringsAsFactors = FALSE)
spx$Date <- as.Date(substr(spx$DateTime, 1, 10))
spx$Close <- as.numeric(spx$Close)
spx <- spx[!is.na(spx$Date) & !is.na(spx$Close), ]

cv_main <- load_main_critical_values(file.path(release_root, "simulations", "main"))

all_results <- list()
for (m_value in c(75L, 100L)) {
  fd_built <- build_fd_data_for_training_window(spx, m = m_value, max_T = 5L)
  all_results[[as.character(m_value)]] <- evaluate_training_window(
    fd_built$fd_data,
    fd_built$dates,
    cv_main,
    m = m_value
  )
}

results <- do.call(rbind, all_results)
summary_table <- build_summary_table(results)

out_dir <- file.path(empirical_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out_dir, "empirical_training_window_robustness_results.csv"), row.names = FALSE, na = "")
write.csv(summary_table, file.path(out_dir, "empirical_training_window_robustness_table.csv"), row.names = FALSE, na = "")
write_latex_table(summary_table, file.path(out_dir, "empirical_training_window_robustness_table.tex"))

message("Saved training-window robustness diagnostics to: ", out_dir)
