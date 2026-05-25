args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
root <- if (length(file_arg) > 0L) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L]))
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

r_files <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("(^|/)outputs?/", r_files)]

failures <- character()
for (path in r_files) {
  ok <- tryCatch({
    parse(path, keep.source = FALSE)
    TRUE
  }, error = function(e) {
    failures <<- c(failures, sprintf("%s :: %s", path, conditionMessage(e)))
    FALSE
  })
  invisible(ok)
}

message("Parsed R files: ", length(r_files))
if (length(failures) > 0L) {
  message(paste(failures, collapse = "\n"))
  quit(status = 1L)
}
