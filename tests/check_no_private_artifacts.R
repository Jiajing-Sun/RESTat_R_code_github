args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
root <- if (length(file_arg) > 0L) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L]))
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

all_files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
rel <- sub(paste0("^", gsub("([].[^$*+?{}|()\\\\])", "\\\\\\1", root), "/?"), "", all_files)
rel <- rel[!grepl("^(\\.git|codex_reports|codex_prompts)(/|$)", rel)]

bad_patterns <- c(
  "^\\.DS_Store$",
  "/\\.DS_Store$",
  "^\\._",
  "/\\._",
  "^__MACOSX/",
  "\\.zip$",
  "\\.RData$",
  "\\.rda$",
  "\\.rds$",
  "\\.RDS$",
  "(^|/)SPX\\.csv$",
  "(^|/)TickData",
  "(^|/)archive/"
)

bad <- unique(rel[vapply(rel, function(x) any(grepl(paste(bad_patterns, collapse = "|"), x)), logical(1L))])
if (length(bad) > 0L) {
  message("Private/generated artifacts found:\n", paste(bad, collapse = "\n"))
  quit(status = 1L)
}

message("No private/generated artifacts detected.")
