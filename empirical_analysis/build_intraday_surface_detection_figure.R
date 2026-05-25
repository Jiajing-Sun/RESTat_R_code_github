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

facet_colours_from_surface <- function(z, palette_fun, alpha = 0.94) {
  z_mid <- 0.25 * (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
  z_lim <- as.numeric(stats::quantile(abs(z_mid), probs = 0.96, na.rm = TRUE))
  if (!is.finite(z_lim) || z_lim <= 0) {
    z_lim <- 1
  }
  scaled <- (z_mid + z_lim) / (2 * z_lim)
  scaled <- pmin(1, pmax(0, scaled))
  idx <- 1 + floor(scaled * 199)
  cols <- palette_fun(200)[idx]
  list(
    cols = grDevices::adjustcolor(cols, alpha.f = alpha),
    z_lim = z_lim
  )
}

draw_detection_plane <- function(pmat, x_range, y0, z_range, fill_col, border_col = fill_col, alpha = 0.15, lty = 1, lwd = 2.2) {
  plane <- trans3d(
    c(x_range[1], x_range[2], x_range[2], x_range[1]),
    rep(y0, 4),
    c(z_range[1], z_range[1], z_range[2], z_range[2]),
    pmat
  )
  polygon(plane$x, plane$y, col = grDevices::adjustcolor(fill_col, alpha.f = alpha), border = NA)

  bottom_edge <- trans3d(x_range, c(y0, y0), c(z_range[1], z_range[1]), pmat)
  lines(bottom_edge$x, bottom_edge$y, col = border_col, lty = lty, lwd = lwd)

  top_edge <- trans3d(x_range, c(y0, y0), c(z_range[2], z_range[2]), pmat)
  lines(top_edge$x, top_edge$y, col = grDevices::adjustcolor(border_col, alpha.f = 0.65), lty = lty, lwd = 1.2)
}

draw_floor_marker <- function(pmat, x_range, y0, z0, col, label = NULL, lty = 1, lwd = 2.8, cex = 0.92, halo_col = grDevices::adjustcolor("white", 0.82)) {
  seg <- trans3d(x_range, c(y0, y0), c(z0, z0), pmat)
  lines(seg$x, seg$y, col = halo_col, lty = lty, lwd = lwd + 2.0)
  lines(seg$x, seg$y, col = col, lty = lty, lwd = lwd)
  if (!is.null(label) && nzchar(label)) {
    lab <- trans3d(max(x_range) + 0.08, y0, z0, pmat)
    text(lab$x, lab$y, labels = label, pos = 4, cex = cex, col = col, xpd = NA)
  }
}

draw_surface_color_key <- function(z_lim, palette_fun) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(fig = c(0.84, 0.96, 0.51, 0.84), new = TRUE, mar = c(1.5, 0.4, 1.2, 3.4), bg = NA)
  z_seq <- seq(-z_lim, z_lim, length.out = 200)
  image(
    x = 1,
    y = z_seq,
    z = matrix(z_seq, nrow = 1),
    col = palette_fun(200),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(4, at = pretty(c(-z_lim, z_lim), n = 5), las = 1, cex.axis = 0.85, col.axis = "gray20", col = "gray55")
  box(col = "gray55")
  mtext("Color scale", side = 3, line = 0.55, cex = 0.88, font = 2, col = "gray20")
}

draw_trace_legend <- function(marker_info) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(fig = c(0.12, 0.88, 0.005, 0.14), new = TRUE, mar = c(0, 0, 0, 0), bg = NA)
  plot.new()
  legend(
    "center",
    legend = marker_info$legend_label,
    col = marker_info$col,
    lty = marker_info$lty,
    lwd = 3.6,
    cex = 0.80,
    bty = "o",
    box.col = "gray70",
    pt.cex = 1.2,
    y.intersp = 1.05,
    text.col = "gray15",
    bg = grDevices::adjustcolor("white", alpha.f = 0.82),
    x.intersp = 0.8,
    ncol = 2
  )
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

plot_start <- as.Date("2020-01-02")
plot_end <- as.Date("2020-04-02")
spx_sub <- spx[spx$Date >= plot_start & spx$Date <= plot_end, , drop = FALSE]

grid_out <- seq(0, 1, length.out = 181)
daily_curves <- lapply(split(spx_sub$Close, spx_sub$Date), smooth_day_curve, grid_out = grid_out)
curve_dates <- as.Date(names(daily_curves))
curve_mat <- do.call(cbind, daily_curves)
z_bp <- 10000 * curve_mat
x_vals <- grid_out
y_vals <- seq_along(curve_dates)

palette_fun <- function(n) grDevices::colorRampPalette(c("#163A70", "#3C8DC4", "#F6F6F1", "#F2A65A", "#B22222"))(n)
surface_cols <- facet_colours_from_surface(z_bp, palette_fun)
facet_cols <- surface_cols$cols
z_color_lim <- surface_cols$z_lim
z_floor <- min(z_bp, na.rm = TRUE) - 2.5
z_top <- max(z_bp, na.rm = TRUE) + 1.0

marker_info <- data.frame(
  legend_label = c(
    "RSMS multiscale MOSUM: Jan 31",
    "RSMS weighted CUSUM: Feb 24",
    "RSMS KS: Feb 27",
    "Stress landmark: Mar 16",
    "HAC KS: Apr 1"
  ),
  date = as.Date(c("2020-01-31", "2020-02-24", "2020-02-27", "2020-03-16", "2020-04-01")),
  col = c("#7B3294", "#D95F02", "#1F78B4", "#202020", "#1B9E77"),
  lty = c(1, 1, 1, 2, 1),
  stringsAsFactors = FALSE
)
marker_info$y <- match(marker_info$date, curve_dates)
marker_info <- marker_info[!is.na(marker_info$y), , drop = FALSE]

output_files <- c(
  file.path(paper_fig_dir, "SPX_2020_IntradaySurface_DetectionDates.png"),
  file.path(bundle_fig_dir, "SPX_2020_IntradaySurface_DetectionDates.png")
)

for (outfile in output_files) {
  png(filename = outfile, width = 2800, height = 2050, res = 240, bg = "#FBFCFF")
  par(fig = c(0, 1, 0.14, 1), mar = c(0.6, 0.35, 3.9, 0.35), xpd = NA, bg = "#FBFCFF")
  pmat <- persp(
    x = x_vals,
    y = y_vals,
    z = z_bp,
    theta = 42,
    phi = 26,
    r = 7.5,
    d = 6,
    expand = 0.78,
    shade = 0.22,
    ltheta = 135,
    lphi = 28,
    col = facet_cols,
    border = grDevices::adjustcolor("#F7F7F7", alpha.f = 0.18),
    ticktype = "detailed",
    nticks = 5,
    xlab = "Standardized Intraday Time",
    ylab = "Trading Day in 2020",
    zlab = "Smoothed 1-minute Return (bp)",
    zlim = c(z_floor, z_top)
  )

  title(main = "Intraday Return Surface with Selected Detection Dates", cex.main = 1.38, font.main = 2, line = 1.05, col.main = "#111111")

  # Mark selected detection dates with semi-transparent planes and floor markers.
  for (ii in order(marker_info$y, decreasing = TRUE)) {
    y0 <- marker_info$y[ii]
    draw_detection_plane(
      pmat,
      x_range = c(min(x_vals), max(x_vals)),
      y0 = y0,
      z_range = c(z_floor, z_top),
      fill_col = marker_info$col[ii],
      border_col = marker_info$col[ii],
      alpha = ifelse(marker_info$lty[ii] == 2, 0.10, 0.14),
      lty = marker_info$lty[ii],
      lwd = 2.2
    )
    draw_floor_marker(
      pmat,
      x_range = c(min(x_vals), max(x_vals)),
      y0 = y0,
      z0 = z_floor,
      col = marker_info$col[ii],
      label = "",
      lty = marker_info$lty[ii],
      lwd = 2.8
    )
  }

  # Add a few explicit date labels along the time axis for orientation.
  axis_dates <- as.Date(c("2020-01-02", "2020-01-31", "2020-02-24", "2020-02-27", "2020-03-16", "2020-04-01"))
  axis_pos <- match(axis_dates, curve_dates)
  axis_lab <- c("Jan 2", "Jan 31", "Feb 24", "Feb 27", "Mar 16", "Apr 1")
  axis_x <- c(-0.06, -0.06, -0.09, -0.02, -0.06, -0.06)
  axis_y_shift <- c(0, 0, -2.2, 2.4, 0, 0)
  axis_text_pos <- c(2, 2, 2, 4, 2, 2)
  axis_cex <- c(0.90, 0.90, 0.88, 0.88, 0.90, 0.90)
  keep <- is.finite(axis_pos)
  for (ii in which(keep)) {
    lab <- trans3d(axis_x[ii], axis_pos[ii] + axis_y_shift[ii], z_floor, pmat)
    text(
      lab$x, lab$y,
      labels = axis_lab[ii],
      cex = axis_cex[ii],
      pos = axis_text_pos[ii],
      col = "gray20",
      xpd = NA
    )
  }

  draw_surface_color_key(z_color_lim, palette_fun)

  mtext(
    "The colored planes mark selected detector dates; the warmer surface tones reveal stress building before the later HAC alarm.",
    side = 3, line = -0.35, cex = 0.98, col = "#2B2B2B"
  )
  draw_trace_legend(marker_info)
  dev.off()
}

message("Saved 3D intraday surface figure to: ", output_files[1])
