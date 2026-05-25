# ==============================================================
# utils.R -- helpers for simulation, parallelism, labels, numerics,
# and Windows-safe atomic file writing
# ============================================================== 

is_windows <- function() identical(.Platform$OS.type, "windows")

safe_solve <- function(A, ridge = 1e-10) {
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  q <- nrow(A)
  A2 <- (A + t(A)) / 2 + ridge * diag(q)
  tryCatch(chol2inv(chol(A2)), error = function(e) solve(A2))
}

safe_matrix <- function(x, ncol = NULL) {
  if (is.null(dim(x))) {
    if (is.null(ncol)) ncol <- 1L
    matrix(x, ncol = ncol)
  } else {
    as.matrix(x)
  }
}

get_connection_limit <- function() {
  opt <- getOption("connections")
  if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt >= 16) return(as.integer(opt))
  envv <- suppressWarnings(as.integer(Sys.getenv("R_CONNECTIONS_LIMIT", unset = "")))
  if (!is.na(envv) && envv >= 16L) return(envv)
  128L
}

current_open_connections <- function() {
  out <- tryCatch(showConnections(all = TRUE), error = function(e) NULL)
  if (is.null(out)) return(0L)
  as.integer(nrow(out))
}

resolve_safe_ncores <- function(requested = NULL,
                                reserve_connections = 10L,
                                hard_cap_windows = 12L,
                                hard_cap_other = 32L) {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1L) detected <- 1L
  detected <- as.integer(detected)

  requested_default <- max(1L, min(detected - 1L, if (is_windows()) hard_cap_windows else hard_cap_other))
  requested_value <- if (is.null(requested) || length(requested) == 0L || is.na(requested)) {
    requested_default
  } else {
    max(1L, as.integer(requested[1L]))
  }

  connection_limit <- get_connection_limit()
  open_connections <- current_open_connections()
  connection_budget <- max(1L, connection_limit - open_connections - as.integer(reserve_connections))
  os_cap <- if (is_windows()) as.integer(hard_cap_windows) else as.integer(hard_cap_other)
  used <- max(1L, min(requested_value, detected, connection_budget, os_cap))

  list(
    requested = as.integer(requested_value),
    detected = as.integer(detected),
    used = as.integer(used),
    connection_limit = as.integer(connection_limit),
    open_connections = as.integer(open_connections),
    reserve_connections = as.integer(reserve_connections),
    connection_budget = as.integer(connection_budget),
    os_cap = as.integer(os_cap)
  )
}

make_psock_cluster <- function(ncores) {
  ncores <- as.integer(max(1L, ncores))
  if (ncores <= 1L) return(NULL)
  parallel::makeCluster(ncores, type = "PSOCK", outfile = "")
}

stop_psock_cluster <- function(cl) {
  if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  invisible(NULL)
}

split_into_chunks <- function(x, nchunks) {
  n <- length(x)
  if (n == 0L) return(list())
  nchunks <- max(1L, min(as.integer(nchunks), n))
  split(x, cut(seq_len(n), breaks = nchunks, labels = FALSE))
}

make_seed_stream <- function(nrep, seed = 123456789L) {
  max_int <- .Machine$integer.max
  seeds <- (as.integer(seed) + seq_len(nrep)) %% max_int
  seeds[seeds == 0L] <- 1L
  seeds
}

canonical_T_scalar <- function(T) {
  if (length(T) != 1L) stop("T must be scalar.")
  if (is.character(T)) {
    x <- trimws(T)
    if (tolower(x) %in% c("inf", "infinity")) return("Inf")
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)) T <- num else return(x)
  }
  if (is.infinite(T)) return("Inf")
  if (!is.numeric(T) || is.na(T)) stop("Invalid T value.")
  if (abs(T - round(T)) < 1e-10) return(as.character(as.integer(round(T))))
  as.character(T)
}

canonical_T_vec <- function(T) vapply(T, canonical_T_scalar, character(1L))

format_q_grid <- function(q_grid) {
  q_grid <- sort(unique(as.integer(q_grid)))
  if (length(q_grid) == 0L) return("<empty>")
  if (length(q_grid) <= 10L) return(paste(q_grid, collapse = ", "))
  sprintf("%d:%d (%d values)", min(q_grid), max(q_grid), length(q_grid))
}

sanitize_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

first_crossing <- function(path, critical_value, k_index = NULL) {
  idx <- which(path > critical_value)
  if (length(idx) == 0L) return(NA_integer_)
  if (is.null(k_index)) return(as.integer(idx[1L]))
  as.integer(k_index[idx[1L]])
}

g_gamma_sq <- function(s, gamma) (1 + s)^2 * (s / (1 + s))^(2 * gamma)

trapz_nonuniform <- function(x, y) {
  n <- length(x)
  if (length(y) != n || n < 2L) return(0)
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

hac_bandwidth_default <- function(m) {
  max(1L, floor(4 * (m / 100)^(2 / 9)))
}

hac_bartlett <- function(Z, bandwidth = NULL, center = TRUE, ridge = 1e-10) {
  Z <- safe_matrix(Z)
  m <- nrow(Z)
  if (m < 2L) return(diag(rep(1, ncol(Z))))
  if (center) Z <- sweep(Z, 2, colMeans(Z), FUN = "-")
  L <- if (is.null(bandwidth)) hac_bandwidth_default(m) else as.integer(max(1L, bandwidth[1L]))
  L <- min(L, m - 1L)
  S <- crossprod(Z) / m
  if (L >= 1L) {
    for (ell in seq_len(L)) {
      w <- 1 - ell / (L + 1)
      Gell <- crossprod(Z[(ell + 1L):m, , drop = FALSE], Z[1L:(m - ell), , drop = FALSE]) / m
      S <- S + w * (Gell + t(Gell))
    }
  }
  S <- (S + t(S)) / 2 + ridge * diag(ncol(Z))
  S
}

build_linear_index_grid <- function(n, size) {
  if (n <= 0L) return(integer(0L))
  size <- max(1L, min(as.integer(size), n))
  unique(as.integer(round(seq(1, n, length.out = size))))
}

build_geometric_index_grid <- function(n, size) {
  if (n <= 0L) return(integer(0L))
  size <- max(1L, min(as.integer(size), n))
  if (size == 1L) return(n)
  g <- exp(seq(log(1), log(n), length.out = size))
  unique(as.integer(round(g)))
}

build_eval_index_grid <- function(k_max, size) {
  out <- build_linear_index_grid(k_max, size)
  out[out >= 1L & out <= k_max]
}

full_or_geometric_lag_grid <- function(k_max, size, exact = FALSE) {
  if (isTRUE(exact) || k_max <= size) return(seq_len(k_max))
  out <- build_geometric_index_grid(k_max, size)
  out[out >= 1L & out <= k_max]
}

normalize_omega_name <- function(x) {
  w <- toupper(trimws(as.character(x)[1L]))
  if (w %in% c("INVSQRT", "INV_SQRT", "1/SQRT", "ONEOVERROOT", "ONEOVERSQRT")) return("InvSqrt")
  stop("Unknown weighted-CUSUM omega. Supported: InvSqrt")
}

make_weighted_cusum_omega <- function(lens, name = "InvSqrt") {
  nm <- normalize_omega_name(name)
  lens <- pmax(as.numeric(lens), 1e-8)
  if (nm == "InvSqrt") return(lens^(-1/2))
  stop("Unknown omega.")
}

normalize_scale_weight_name <- function(x) {
  w <- toupper(trimws(as.character(x)[1L]))
  if (w %in% c("EQUAL", "UNIFORM", "ONE")) return("Equal")
  stop("Unknown multiscale weight. Supported: Equal")
}

make_multiscale_weight <- function(h, name = "Equal") {
  nm <- normalize_scale_weight_name(name)
  if (nm == "Equal") return(rep(1, length(h)))
  stop("Unknown multiscale weight.")
}

normalize_weight_name <- function(weight) {
  w <- toupper(trimws(as.character(weight)[1L]))
  if (w %in% c("U", "UNIFORM", "CONST", "CONSTANT", "ONE", "WU", "W_U", "1")) return("U")
  if (w %in% c("EARLY", "W_EARLY")) return("Early")
  if (w %in% c("MID", "W_MID")) return("Mid")
  if (w %in% c("LATE", "W_LATE")) return("Late")
  stop("Unknown CvM weight. Use one of: U, Early, Mid, Late.")
}

make_cvm_weight_finite <- function(k_vec, m, T, weight = "U") {
  w <- normalize_weight_name(weight)
  tau <- (k_vec / m) / T
  if (w == "U") return(rep(1, length(k_vec)))
  if (w == "Late") return(2 * tau)
  if (w == "Early") return(2 * (1 - tau))
  if (w == "Mid") return(6 * tau * (1 - tau))
  stop("Unknown finite-horizon weight.")
}

leading_mahalanobis_sequence <- function(D, V, ridge = 1e-10) {
  D <- safe_matrix(D)
  q_max <- ncol(D)
  out <- matrix(NA_real_, nrow = nrow(D), ncol = q_max)
  for (q in seq_len(q_max)) {
    Vinv <- safe_solve(V[1:q, 1:q, drop = FALSE], ridge = ridge)
    Dq <- D[, 1:q, drop = FALSE]
    out[, q] <- rowSums((Dq %*% Vinv) * Dq)
  }
  out
}

row_cumsum_weighted_squares <- function(D, w = NULL) {
  D <- safe_matrix(D)
  if (is.null(w)) w <- rep(1, ncol(D))
  w <- as.numeric(w)
  tmp <- sweep(D^2, 2, w, FUN = "*")
  t(apply(tmp, 1, cumsum))
}

method_label_main <- function(standardizer, type, gamma = NA_real_, weight_name = NA_character_) {
  if (toupper(type) == "KS") {
    sprintf("%s KS (gamma=%s)", standardizer, format(gamma, trim = TRUE))
  } else {
    sprintf("%s weighted CvM [%s]", standardizer, normalize_weight_name(weight_name))
  }
}

method_label_alt <- function(standardizer, detector, gamma = NA_real_, bandwidth_h = NA_real_, omega_name = NA_character_, hset_name = NA_character_, scale_weight_name = NA_character_) {
  if (detector %in% c("FullCUSUM", "PageCUSUM")) {
    sprintf("%s %s (gamma=%s)", standardizer, detector, format(gamma, trim = TRUE))
  } else if (detector == "WeightedCUSUM") {
    sprintf("%s %s (gamma=%s, %s)", standardizer, detector, format(gamma, trim = TRUE), omega_name)
  } else if (detector == "MOSUM") {
    sprintf("%s MOSUM (h=%s)", standardizer, format(bandwidth_h, trim = TRUE))
  } else if (detector == "MultiscaleMOSUM") {
    sprintf("%s MultiscaleMOSUM (%s, %s)", standardizer, hset_name, scale_weight_name)
  } else {
    sprintf("%s %s", standardizer, detector)
  }
}

method_group_from_row <- function(family, type, detector) {
  if (family == "Main" && toupper(type) == "KS") return("Main_KS")
  if (family == "Main" && toupper(type) == "CVM") return("Main_CvM")
  if (detector %in% c("FullCUSUM", "PageCUSUM", "WeightedCUSUM")) return("Benchmark_CUSUM")
  if (detector %in% c("MOSUM", "MultiscaleMOSUM")) return("Benchmark_MOSUM")
  "Other"
}

file_nonempty <- function(path, min_bytes = 8L) {
  if (length(path) == 0L) return(logical(0))
  exists <- file.exists(path)
  size <- suppressWarnings(file.info(path)$size)
  exists & is.finite(size) & (size >= min_bytes)
}

safe_unlink <- function(path) {
  if (length(path) == 0L || is.na(path)) return(invisible(FALSE))
  try(unlink(path, force = TRUE), silent = TRUE)
  invisible(TRUE)
}

atomic_save_rds <- function(obj, path, compress = "gzip") {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "_tmp_"), tmpdir = dirname(path), fileext = ".rds")
  on.exit(safe_unlink(tmp), add = TRUE)
  saveRDS(obj, file = tmp, compress = compress)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    safe_unlink(tmp)
  }
  if (!ok) stop("Failed to atomically write RDS to: ", path)
  invisible(path)
}

atomic_write_csv <- function(df, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "_tmp_"), tmpdir = dirname(path), fileext = ".csv")
  on.exit(safe_unlink(tmp), add = TRUE)
  utils::write.csv(df, file = tmp, row.names = FALSE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    safe_unlink(tmp)
  }
  if (!ok) stop("Failed to atomically write CSV to: ", path)
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "_tmp_"), tmpdir = dirname(path), fileext = ".txt")
  on.exit(safe_unlink(tmp), add = TRUE)
  writeLines(lines, con = tmp, useBytes = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    safe_unlink(tmp)
  }
  if (!ok) stop("Failed to atomically write text file to: ", path)
  invisible(path)
}

safe_read_rds <- function(path) {
  if (!file_nonempty(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

safe_read_csv <- function(path) {
  if (!file_nonempty(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

combine_rbind <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1L))]
  if (length(lst) == 0L) return(data.frame())
  do.call(rbind, lst)
}

list_batch_rds <- function(batch_dir) {
  if (!dir.exists(batch_dir)) return(character(0))
  sort(list.files(batch_dir, pattern = "^batch_[0-9]+\\.rds$", full.names = TRUE))
}

append_progress_line <- function(path, text) {
  ensure_dir(dirname(path))
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text, "\n"),
      file = path, append = TRUE)
  invisible(path)
}

write_run_summary <- function(lines, output_dir, filename = "run_summary.txt") {
  ensure_dir(output_dir)
  atomic_write_lines(lines, file.path(output_dir, filename))
}
