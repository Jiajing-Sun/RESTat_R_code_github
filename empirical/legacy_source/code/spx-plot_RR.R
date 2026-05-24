library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

# --- Set working directory (自动化) ---
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("--file=", args, value = TRUE)
  if (length(file.arg) > 0) {
    script.path <- normalizePath(sub("--file=", "", file.arg))
    setwd(dirname(script.path))
  }
}

# Load data
spx <- read_csv("SPX.csv", show_col_types = FALSE) %>%
  mutate(
    Date = as_date(DateTime),
    Year = year(Date)
  ) %>%
  filter(Year == 2020) %>%
  arrange(Date, DateTime)

# 1-min log returns within each day
spx_ret <- spx %>%
  group_by(Date) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(r = log(Close) - log(lag(Close))) %>%
  filter(!is.na(r)) %>%
  ungroup()

# Daily return: sum of intraday log returns
spx_daily_ret <- spx_ret %>%
  group_by(Date) %>%
  summarise(DailyLogRet = sum(r, na.rm = TRUE), .groups = "drop")

# COVID crash date
crash_date <- ymd("2020-03-16")

# Plot
p <- ggplot(spx_daily_ret, aes(x = Date, y = DailyLogRet)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(
    xintercept = as.numeric(crash_date),
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  labs(
    x = "Date",
    y = "Daily log return",
    title = "S&P 500 Daily Log Returns in 2020"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )


ggsave("SPX_2020_DailyLogReturn.png", p, width = 6, height = 4, units = "in", dpi = 1200)
