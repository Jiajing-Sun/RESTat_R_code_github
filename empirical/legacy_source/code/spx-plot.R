# Load libraries
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(ggplot2)

# Set folder
# --- Set working directory (自动化) ---
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  script.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script.dir)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("--file=", args, value = TRUE)
  if (length(file.arg) > 0) {
    script.path <- normalizePath(sub("--file=", "", file.arg))
    setwd(dirname(script.path))
  }
}

# Load data
spx_data <- read_csv("SPX.csv", show_col_types = FALSE)

head(spx_data)


# Preprocess: extract 2020 and compute daily average
spx_2020 <- spx_data %>%
  mutate(Date = as_date(DateTime)) %>%
  filter(year(Date) == 2020) %>%
  group_by(Date) %>%
  summarise(DailyAvg = mean(Close, na.rm = TRUE), .groups = "drop")

# COVID-19 market crash date
crash_date <- ymd("2020-03-16")

# Plot
p <- ggplot(spx_2020, aes(x = Date, y = DailyAvg)) +
  geom_line(color = "black", linewidth = 0.4) +
  geom_vline(xintercept = as.numeric(crash_date), linetype = "dashed", color = "black") +
  labs(
    x = "Date",
    y = "Average Close Price",
    title = "S&P 500 Daily Average Closing Price in 2020"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

ggsave("SPX_2020_COVID_Crash.png", p, width = 6, height = 4, units = "in", dpi = 1200)
