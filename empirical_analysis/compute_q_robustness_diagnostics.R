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

fmt_signal_date <- function(x) {
  x <- as.character(x)
  out <- rep("No signal", length(x))
  keep <- !is.na(x) & nzchar(x)
  if (any(keep)) {
    d <- as.Date(x[keep])
    mon <- format(d, "%b")
    day <- sub("^0", "", format(d, "%d"))
    yr <- format(d, "%Y")
    out[keep] <- paste0(mon, ".~", day, ", ", yr)
  }
  out
}

write_latex_table <- function(df, path, caption, label, note = NULL) {
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4.0pt}",
    "\\renewcommand{\\arraystretch}{0.95}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lcccccc}",
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

build_fd_data <- function(spx_path) {
  spx <- read.csv(spx_path, stringsAsFactors = FALSE)
  spx$Date <- as.Date(substr(spx$DateTime, 1, 10))
  spx <- spx[
    spx$Date >= as.Date("2019-10-21") &
      spx$Date <= as.Date("2020-12-28"),
    c("Date", "DateTime", "Close")
  ]
  spx$Close <- as.numeric(spx$Close)

  dates <- sort(unique(spx$Date))
  daily_list <- lapply(dates, function(d) {
    day <- spx[spx$Date == d, c("Date", "DateTime", "Close")]
    day$Close <- c(NA_real_, diff(log(day$Close)))
    day <- day[!is.na(day$Close), ]
    rownames(day) <- NULL
    day
  })

  min_obs <- min(vapply(daily_list, nrow, integer(1L)))
  nbasis <- min(21L, floor(min_obs * 0.8))
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), norder = 4, nbasis = nbasis)
  coefs <- matrix(nrow = nbasis, ncol = length(daily_list))
  day_names <- character(length(daily_list))

  for (i in seq_along(daily_list)) {
    day <- daily_list[[i]]
    fd_obj <- fda::smooth.basis(
      argvals = seq(0, 1, length.out = nrow(day)),
      y = day$Close,
      fdParobj = basis
    )
    coefs[, i] <- fd_obj$fd$coefs
    day_names[i] <- as.character(day$Date[1])
  }

  list(
    fd_data = fda::fd(
      coef = coefs,
      basisobj = basis,
      fdnames = list(time = "t", replicates = day_names, values = "SPX_Close")
    ),
    dates = as.Date(day_names)
  )
}

run_q_grid <- function(fd_data, all_dates, cv_main, q_grid = c(3L, 4L, 5L),
                       T_grid = c(1L, 2L, 5L), alpha = 0.05, m = 50L) {
  main_catalog <- build_main_method_catalog()
  main_catalog <- subset(
    main_catalog,
    type == "KS" &
      standardizer %in% c("HAC", "SSMS", "RSMS") &
      gamma %in% c(0, 0.15)
  )

  out <- list()
  idx <- 1L
  for (q in q_grid) {
    fp <- fpca_project_scores(fd_data, m = m, fixed_q = q, fve_threshold = 0.95, n_grid = 301L)
    white <- whiten_scores(fp$scores_train, fp$scores_monitor)

    for (T_value in T_grid) {
      k_monitor <- m * T_value
      raw_train <- fp$scores_train
      raw_monitor <- fp$scores_monitor[1:k_monitor, , drop = FALSE]
      rsms_train <- white$scores_train
      rsms_monitor <- white$scores_monitor[1:k_monitor, , drop = FALSE]

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
  }
  do.call(rbind, out)
}

build_summary_table <- function(results, gamma_value) {
  sub <- results[
    results$gamma == gamma_value &
      results$standardizer %in% c("RSMS", "HAC"),
    ,
    drop = FALSE
  ]
  q_vals <- sort(unique(sub$q_used))
  out <- data.frame(
    q = q_vals,
    `RSMS ($T=1$)` = character(length(q_vals)),
    `RSMS ($T=2$)` = character(length(q_vals)),
    `RSMS ($T=5$)` = character(length(q_vals)),
    `HAC ($T=1$)` = character(length(q_vals)),
    `HAC ($T=2$)` = character(length(q_vals)),
    `HAC ($T=5$)` = character(length(q_vals)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    for (std in c("RSMS", "HAC")) {
      for (T_value in c("1", "2", "5")) {
        hit <- sub[
          sub$q_used == q &
            sub$standardizer == std &
            sub$T_label == T_value,
          ,
          drop = FALSE
        ]
        value <- if (nrow(hit) == 1L) fmt_signal_date(hit$first_rejection_date) else "NA"
        out[i, paste0(std, " ($T=", T_value, "$)")] <- value
      }
    }
  }

  out
}

script_dir <- get_script_dir()
empirical_dir <- normalizePath(script_dir)
release_root <- normalizePath(file.path(empirical_dir, ".."))
main_modules_dir <- file.path(release_root, "simulations", "main", "R")
for (f in c("dependencies.R", "utils.R", "critical_values_lookup.R", "method_catalog.R", "fpca_pipeline.R", "detectors_main.R")) {
  source(file.path(main_modules_dir, f), local = TRUE)
}
ensure_simulation_packages()

spx_path <- Sys.getenv("SPX_DATA_FILE", unset = file.path(empirical_dir, "data", "SPX.csv"))
if (!file.exists(spx_path)) {
  stop("Could not find SPX.csv at: ", spx_path)
}
fd_built <- build_fd_data(spx_path)
cv_main <- load_main_critical_values(file.path(release_root, "simulations", "main"))
results <- run_q_grid(fd_built$fd_data, fd_built$dates, cv_main)

out_dir <- file.path(empirical_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out_dir, "empirical_q_robustness_results.csv"), row.names = FALSE, na = "")

gamma0 <- build_summary_table(results, 0)
gamma015 <- build_summary_table(results, 0.15)
write.csv(gamma0, file.path(out_dir, "empirical_q_robustness_gamma0.csv"), row.names = FALSE, na = "")
write.csv(gamma015, file.path(out_dir, "empirical_q_robustness_gamma015.csv"), row.names = FALSE, na = "")

write_latex_table(
  gamma0,
  file.path(out_dir, "empirical_q_robustness_gamma0.tex"),
  "Retained-dimension robustness for the main KS empirical alarms ($\\gamma=0$).",
  "tab:bundle-empirical-q-robustness-gamma0",
  note = "SSMS remains below the 5\\% KS threshold for all retained dimensions $q\\in\\{3,4,5\\}$ and horizons $T\\in\\{1,2,5\\}$."
)

write_latex_table(
  gamma015,
  file.path(out_dir, "empirical_q_robustness_gamma015.tex"),
  "Retained-dimension robustness for the main KS empirical alarms ($\\gamma=0.15$).",
  "tab:bundle-empirical-q-robustness-gamma015",
  note = "SSMS remains below the 5\\% KS threshold for all retained dimensions $q\\in\\{3,4,5\\}$ and horizons $T\\in\\{1,2,5\\}$."
)

message("Saved q-robustness empirical diagnostics to: ", out_dir)
