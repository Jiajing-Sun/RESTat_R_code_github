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

smooth_day_curve <- function(day_close, grid_out) {
  if (length(day_close) < 10) {
    return(rep(NA_real_, length(grid_out)))
  }
  ret <- diff(log(day_close))
  t_in <- seq(0, 1, length.out = length(ret))
  fit <- smooth.spline(x = t_in, y = ret, all.knots = FALSE)
  predict(fit, x = grid_out)$y
}

mean_curve_for_window <- function(data, start_date, end_date, grid_out) {
  sub <- data[data$Date >= start_date & data$Date <= end_date, , drop = FALSE]
  if (nrow(sub) == 0) {
    stop("No observations found between ", start_date, " and ", end_date)
  }
  daily_curves <- lapply(split(sub$Close, sub$Date), smooth_day_curve, grid_out = grid_out)
  curve_mat <- do.call(rbind, daily_curves)
  colMeans(curve_mat, na.rm = TRUE)
}

script_dir <- get_script_dir()
empirical_dir <- normalizePath(script_dir)
data_file <- Sys.getenv("SPX_DATA_FILE", unset = file.path(empirical_dir, "data", "SPX.csv"))
bundle_fig_dir <- file.path(empirical_dir, "figures")
paper_fig_dir <- bundle_fig_dir

dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bundle_fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Could not find SPX.csv at: ", data_file)
}
spx <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
blank_names <- names(spx) == ""
if (any(blank_names)) {
  spx <- spx[, !blank_names, drop = FALSE]
}

spx$DateTime <- as.POSIXct(spx$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
spx$Date <- as.Date(spx$DateTime)
spx$Close <- as.numeric(spx$Close)
spx <- spx[!is.na(spx$Date) & !is.na(spx$Close), ]
spx <- spx[order(spx$Date, spx$DateTime), , drop = FALSE]

grid_out <- seq(0, 1, length.out = 301)

training_curve <- mean_curve_for_window(spx, as.Date("2019-10-21"), as.Date("2019-12-31"), grid_out)
transition_curve <- mean_curve_for_window(spx, as.Date("2020-02-24"), as.Date("2020-02-28"), grid_out)
stress_curve <- mean_curve_for_window(spx, as.Date("2020-03-09"), as.Date("2020-03-13"), grid_out)

ylim_all <- range(c(training_curve, transition_curve, stress_curve), finite = TRUE)
ylim_pad <- 0.12 * diff(ylim_all)
ylim_use <- c(ylim_all[1] - ylim_pad, ylim_all[2] + ylim_pad)

regime_colors <- c(
  Training = "#0072B2",
  Transition = "#D55E00",
  Stress = "#CC79A7"
)

output_files <- c(
  file.path(paper_fig_dir, "SPX_2020_IntradayCurve_Regimes.png"),
  file.path(bundle_fig_dir, "SPX_2020_IntradayCurve_Regimes.png")
)

for (outfile in output_files) {
  png(filename = outfile, width = 2000, height = 1300, res = 220)
  par(mar = c(7.8, 8.4, 1.5, 1.2) + 0.1, las = 1, xpd = NA)
  plot(
    grid_out, training_curve,
    type = "l",
    lwd = 4.0,
    lty = 1,
    col = regime_colors[["Training"]],
    ylim = ylim_use,
    xaxt = "n",
    xlab = "",
    ylab = "",
    cex.lab = 1.30,
    cex.axis = 1.1
  )
  lines(grid_out, transition_curve, lwd = 4.0, lty = 1, col = regime_colors[["Transition"]])
  lines(grid_out, stress_curve, lwd = 4.0, lty = 1, col = regime_colors[["Stress"]])
  axis(1, at = c(0, 0.5, 1), labels = c("Open", "Midday", "Close"), cex.axis = 1.15)
  title(xlab = "Standardized Intraday Time", line = 3.2, cex.lab = 1.30)
  title(ylab = "Smoothed 1-minute Return", line = 5.3, cex.lab = 1.30)
  abline(h = 0, lty = 3, lwd = 1.2, col = "gray40")
  legend(
    "bottom",
    inset = c(0, -0.33),
    legend = c("Training period", "Late-February transition", "March stress"),
    col = unname(regime_colors[c("Training", "Transition", "Stress")]),
    lwd = 4.0,
    lty = 1,
    horiz = TRUE,
    bty = "n",
    cex = 0.92,
    x.intersp = 0.65,
    seg.len = 2.0
  )
  dev.off()
}

message("Saved intraday regime figure to: ", output_files[1])
