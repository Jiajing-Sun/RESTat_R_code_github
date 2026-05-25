# Compute the finite-sample off-diagonal diagnostic used to assess the
# approximate diagonalization condition for the RSMS standardizer.

decode_rscript_path <- function(x) gsub("~\\+~", " ", x)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(decode_rscript_path(sub("^--file=", "", file_arg[1L])))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

parse_named_args_local <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (arg in raw) {
    if (!grepl("^--", arg)) next
    key_val <- sub("^--", "", arg)
    parts <- strsplit(key_val, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- if (length(parts) > 1L) paste(parts[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

offdiag_ratio <- function(scores_train, ridge = 1e-10, hac_bandwidth = NULL) {
  scores_train <- safe_matrix(scores_train)
  whitened <- whiten_scores(scores_train, scores_train, ridge = ridge)
  centered <- sweep(scores_train, 2, colMeans(scores_train), FUN = "-")
  gamma_hat <- hac_bartlett(centered, bandwidth = hac_bandwidth, center = FALSE, ridge = ridge)
  w_hat <- whitened$transform
  transformed <- t(w_hat) %*% gamma_hat %*% w_hat
  off <- transformed
  diag(off) <- 0
  denom <- sqrt(sum(transformed^2))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sqrt(sum(off^2)) / denom
}

build_empirical_curve_values <- function(spx_path) {
  spx <- read.csv(spx_path, stringsAsFactors = FALSE, check.names = FALSE)
  blank_names <- names(spx) == ""
  if (any(blank_names)) spx <- spx[, !blank_names, drop = FALSE]
  if (!all(c("DateTime", "Close") %in% names(spx))) {
    stop("SPX.csv must contain DateTime and Close columns.")
  }
  spx$Date <- as.Date(substr(spx$DateTime, 1, 10))
  spx$Close <- as.numeric(spx$Close)
  spx <- spx[
    !is.na(spx$Date) & !is.na(spx$Close) &
      spx$Date >= as.Date("2019-10-21") &
      spx$Date <= as.Date("2020-12-28"),
    c("Date", "DateTime", "Close")
  ]

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
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), norder = 4L, nbasis = nbasis)
  eval_grid <- seq(0, 1, length.out = 301L)
  values <- matrix(NA_real_, nrow = length(daily_list), ncol = length(eval_grid))
  day_names <- character(length(daily_list))

  for (i in seq_along(daily_list)) {
    day <- daily_list[[i]]
    fd_obj <- fda::smooth.basis(
      argvals = seq(0, 1, length.out = nrow(day)),
      y = day$Close,
      fdParobj = basis
    )
    values[i, ] <- as.numeric(fda::eval.fd(eval_grid, fd_obj$fd))
    day_names[i] <- as.character(day$Date[1])
  }

  list(
    curve_values = values,
    dates = as.Date(day_names),
    nbasis = nbasis,
    min_obs = min_obs
  )
}

script_dir <- get_script_dir()
simulations_dir <- normalizePath(file.path(script_dir, ".."))
release_root <- normalizePath(file.path(simulations_dir, ".."))
module_dir <- file.path(simulations_dir, "contaminated_training", "R")

for (f in c(
  "utils.R", "fpca_pipeline.R", "genData.R", "scenarios.R"
)) {
  source(file.path(module_dir, f), local = TRUE)
}

args <- parse_named_args_local()
nsim <- as.integer(args$nsim %||% "1000")
if (is.na(nsim) || nsim < 1L) stop("--nsim must be a positive integer.")
seed0 <- as.integer(args$seed %||% "20260515")
if (is.na(seed0)) seed0 <- 20260515L

out_dir <- file.path(script_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Representative simulation design: fMA(1), m=500, T=2, clean training sample,
# 95 percent FVE selection, 301-point evaluation grid, and default Bartlett HAC.
# The diagnostic is computed on the matrix/SVD implementation used in the
# finite-sample simulation code; it does not require critical-value simulation.
sim_rows <- vector("list", nsim)
for (r in seq_len(nsim)) {
  set.seed(seed0 + r)
  base <- generate_base_curves("fMA1", total_n = 500L + 500L * 2L, n_grid = 301L)
  fp <- fpca_project_scores(
    base$values,
    m = 500L,
    q_cap = 30L,
    fve_threshold = 0.95,
    fixed_q = NA_integer_,
    n_grid = 301L
  )
  ratio <- offdiag_ratio(fp$scores_train, ridge = 1e-10, hac_bandwidth = NULL)
  sim_rows[[r]] <- data.frame(
    replication = r,
    dgp = "fMA(1)",
    m = 500L,
    T = 2L,
    q_used = fp$q_used,
    delta_off = ratio,
    stringsAsFactors = FALSE
  )
  if (r %% 50L == 0L) {
    message(sprintf("Completed %d/%d simulation diagnostic replications", r, nsim))
  }
}
sim_raw <- do.call(rbind, sim_rows)
write.csv(sim_raw, file.path(out_dir, "rsms_offdiag_simulation_raw.csv"), row.names = FALSE)

sim_summary <- data.frame(
  design = "fMA(1), m=500, T=2, 95% FVE",
  replications = nsim,
  mean_delta_off = mean(sim_raw$delta_off, na.rm = TRUE),
  median_delta_off = stats::median(sim_raw$delta_off, na.rm = TRUE),
  p90_delta_off = as.numeric(stats::quantile(sim_raw$delta_off, 0.90, na.rm = TRUE)),
  mean_q = mean(sim_raw$q_used, na.rm = TRUE),
  median_q = stats::median(sim_raw$q_used, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(sim_summary, file.path(out_dir, "rsms_offdiag_simulation_summary.csv"), row.names = FALSE)

spx_path <- Sys.getenv(
  "SPX_DATA_FILE",
  unset = file.path(release_root, "empirical_analysis", "data", "SPX.csv")
)
emp_rows <- list()
if (file.exists(spx_path)) {
  emp <- build_empirical_curve_values(spx_path)
  for (q in c(3L, 4L, 5L)) {
    fp <- fpca_project_scores(
      emp$curve_values,
      m = 50L,
      q_cap = 30L,
      fve_threshold = 0.95,
      fixed_q = q,
      n_grid = 301L
    )
    emp_rows[[length(emp_rows) + 1L]] <- data.frame(
      q = q,
      q_used = fp$q_used,
      m = 50L,
      start_date = as.character(emp$dates[1L]),
      training_end_date = as.character(emp$dates[50L]),
      monitor_start_date = as.character(emp$dates[51L]),
      nbasis = emp$nbasis,
      min_intraday_returns = emp$min_obs,
      delta_off = offdiag_ratio(fp$scores_train, ridge = 1e-10, hac_bandwidth = NULL),
      stringsAsFactors = FALSE
    )
  }
  emp_out <- do.call(rbind, emp_rows)
  write.csv(emp_out, file.path(out_dir, "rsms_offdiag_empirical.csv"), row.names = FALSE)
}

message("Saved RSMS off-diagonal diagnostics to: ", out_dir)
