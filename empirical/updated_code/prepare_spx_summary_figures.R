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

script_dir <- get_script_dir()
empirical_dir <- normalizePath(file.path(script_dir, ".."))
data_file <- file.path(empirical_dir, "legacy_source", "data", "SPX.csv")
fig_dir <- file.path(empirical_dir, "figures")

if (!file.exists(data_file)) {
  stop("Could not find SPX.csv at: ", data_file)
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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
spx_2020 <- spx[format(spx$Date, "%Y") == "2020", ]

if (nrow(spx_2020) == 0) {
  stop("No 2020 observations found in SPX.csv.")
}

daily_avg <- data.frame(
  Date = as.Date(names(tapply(spx_2020$Close, spx_2020$Date, mean, na.rm = TRUE))),
  Close = as.numeric(tapply(spx_2020$Close, spx_2020$Date, mean, na.rm = TRUE))
)
daily_avg <- daily_avg[order(daily_avg$Date), , drop = FALSE]

daily_logret <- do.call(
  rbind,
  lapply(split(spx_2020, spx_2020$Date), function(day_df) {
    day_df <- day_df[order(day_df$DateTime), ]
    if (nrow(day_df) < 2) {
      return(data.frame(Date = unique(day_df$Date), DailyLogRet = NA_real_))
    }
    r <- diff(log(day_df$Close))
    data.frame(Date = unique(day_df$Date), DailyLogRet = sum(r, na.rm = TRUE))
  })
)
daily_logret <- daily_logret[!is.na(daily_logret$DailyLogRet), ]
daily_logret <- daily_logret[order(daily_logret$Date), ]

crash_date <- as.Date("2020-03-16")

png(
  filename = file.path(fig_dir, "SPX_2020_Min_Close.png"),
  width = 1800,
  height = 1200,
  res = 220
)
par(mar = c(5, 5, 3, 1) + 0.1, las = 1, cex.axis = 0.9, cex.lab = 1.0)
plot(
  daily_avg$Date,
  daily_avg$Close,
  type = "l",
  lwd = 1.4,
  col = "black",
  xlab = "Date",
  ylab = "Average Close Price",
  main = "S&P 500 Daily Average Closing Price in 2020"
)
abline(v = crash_date, lty = 2, lwd = 1.1, col = "black")
dev.off()

png(
  filename = file.path(fig_dir, "SPX_2020_DailyLogReturn.png"),
  width = 1800,
  height = 1200,
  res = 220
)
par(mar = c(5, 5, 3, 1) + 0.1, las = 1, cex.axis = 0.9, cex.lab = 1.0)
plot(
  daily_logret$Date,
  daily_logret$DailyLogRet,
  type = "l",
  lwd = 1.2,
  col = "black",
  xlab = "Date",
  ylab = "Daily log return",
  main = "S&P 500 Daily Log Returns in 2020"
)
abline(h = 0, lwd = 1.0, col = "gray40")
abline(v = crash_date, lty = 2, lwd = 1.1, col = "black")
dev.off()

message("Saved figures to: ", fig_dir)
