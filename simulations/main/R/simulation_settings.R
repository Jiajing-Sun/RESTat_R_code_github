# ==============================================================
# simulation_settings.R -- shared defaults, CLI parsing, and
# manuscript-aligned scenario metadata
# ==============================================================

streaming_curve_basis_dimension <- function() 21L

localized_change_basis_index_default <- function() 5L

streaming_curve_supported_dgp_types <- function() c("BB", "fIID", "fMA1")

streaming_curve_supported_power_scenarios <- function() {
  c("level_shift", "smooth_change", "abrupt_local_change", "sinusoidal_change")
}

streaming_curve_scenario_formulas <- function() {
  list(
    level_shift = "Delta",
    smooth_change = "Delta * (t - t*) / (mT)",
    abrupt_local_change = "Delta * b_k(s)",
    sinusoidal_change = "Delta * sin(pi * (t - t*) / (mT))"
  )
}

streaming_curve_hac_bandwidth_formula <- function() "4 * (m / 100)^(2 / 9)"

parse_csv_tokens <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(character(0L))
  raw_tokens <- unlist(strsplit(as.character(x[1L]), ",", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(raw_tokens)
  tokens[nzchar(tokens)]
}

parse_bool_arg <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(isTRUE(default))
  key <- toupper(trimws(as.character(x[1L])))
  if (key %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (key %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop("Invalid logical value: ", x[1L])
}

parse_character_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  unique(tokens)
}

parse_numeric_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- suppressWarnings(as.numeric(tokens))
  if (any(is.na(vals))) stop("Invalid numeric vector argument: ", as.character(x[1L]))
  unique(vals)
}

parse_T_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- vapply(tokens, function(tok) {
    key <- tolower(trimws(tok))
    if (key %in% c("inf", "infinity")) return(Inf)
    num <- suppressWarnings(as.numeric(tok))
    if (is.na(num)) stop("Invalid T-grid entry: ", tok)
    num
  }, numeric(1L))
  unique(vals)
}

parse_integer_vector_arg <- function(x, default = NULL) {
  tokens <- parse_csv_tokens(x)
  if (length(tokens) == 0L) return(default)
  vals <- unlist(lapply(tokens, function(tok) {
    if (grepl("^-?[0-9]+:-?[0-9]+$", tok)) {
      endpoints <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1L]])
      return(seq(endpoints[1L], endpoints[2L]))
    }
    val <- suppressWarnings(as.numeric(tok))
    if (is.na(val) || abs(val - round(val)) > 1e-10) stop("Expected integer value, got: ", tok)
    as.integer(round(val))
  }), use.names = FALSE)
  unique(as.integer(vals))
}

parse_scalar_integer_arg <- function(x, default) {
  vals <- parse_integer_vector_arg(x, default = default)
  as.integer(vals[1L])
}

parse_scalar_numeric_arg <- function(x, default) {
  vals <- parse_numeric_vector_arg(x, default = default)
  as.numeric(vals[1L])
}

validate_nonempty_vector <- function(x, label) {
  if (length(x) == 0L) stop(label, " must not be empty.")
  x
}

validate_supported_values <- function(x, supported, label) {
  bad <- setdiff(x, supported)
  if (length(bad) > 0L) {
    stop(label, " contains unsupported values: ", paste(bad, collapse = ", "))
  }
  x
}

default_shared_simulation_settings <- function() {
  list(
    dgp_types = streaming_curve_supported_dgp_types(),
    alpha_levels = c(0.10, 0.05, 0.01),
    gamma_vec = c(0, 0.15),
    cvm_weights = c("U", "Early", "Mid", "Late"),
    include_alt_detectors = TRUE,
    mosum_h_vec = c(0.10, 0.20),
    weighted_omega_names = c("InvSqrt"),
    multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
    multiscale_scale_names = c("Equal"),
    nsim = 1000L,
    basis_k = localized_change_basis_index_default(),
    q_cap = 30L,
    fve_threshold = 0.95,
    fixed_q = NA_integer_,
    n_grid = 301L,
    ridge = 1e-10,
    range_floor = 1e-8,
    hac_bandwidth = NULL,
    page_length_grid_size = 40L,
    weighted_length_grid_size = 40L,
    finite_eval_grid_size = 250L,
    exact_page_scan = FALSE,
    exact_weighted_scan = FALSE,
    batch_size = 10L,
    overwrite_existing = FALSE
  )
}

default_null_run_settings <- function() {
  modifyList(default_shared_simulation_settings(), list(
    mode = "null",
    m_vals = c(100L, 200L, 500L, 1000L),
    T_grid = c(1, 2, 5, 10)
  ))
}

power_delta_grid_default <- function(scenario, dgp_type, include_zero = TRUE) {
  scenario <- match.arg(scenario, streaming_curve_supported_power_scenarios())
  dgp_type <- match.arg(dgp_type, streaming_curve_supported_dgp_types())

  base_vals <- if (scenario %in% c("level_shift", "smooth_change", "sinusoidal_change")) {
    if (dgp_type == "BB") {
      c(0.050, 0.058, 0.065, 0.073, 0.080, 0.088, 0.095, 0.103, 0.110, 0.118, 0.125)
    } else {
      c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029)
    }
  } else if (dgp_type == "BB") {
    c(0.100, 0.400, 0.700, 1.000, 1.300, 1.600)
  } else {
    c(0.05, 0.08, 0.11, 0.14, 0.17, 0.20)
  }

  if (isTRUE(include_zero)) c(0, base_vals) else base_vals
}

default_power_run_settings <- function(scenario = "level_shift") {
  scenario <- match.arg(scenario, streaming_curve_supported_power_scenarios())
  shared <- default_shared_simulation_settings()
  dgp_types <- streaming_curve_supported_dgp_types()
  delta_map <- stats::setNames(vector("list", length(dgp_types)), dgp_types)
  for (nm in dgp_types) delta_map[[nm]] <- power_delta_grid_default(scenario, nm, include_zero = TRUE)

  modifyList(shared, list(
    mode = "power",
    scenario = scenario,
    m_vals = 500L,
    T_grid = c(1, 2, 5, 10),
    s_star_vals = c(50L, 200L),
    delta_map = delta_map
  ))
}

format_delta_map_summary <- function(delta_map) {
  parts <- vapply(names(delta_map), function(nm) {
    vals <- delta_map[[nm]]
    sprintf("%s: %s", nm, paste(format(vals, trim = TRUE, scientific = FALSE), collapse = ", "))
  }, character(1L))
  paste(parts, collapse = " | ")
}

resolve_null_run_settings <- function(args) {
  defaults <- default_null_run_settings()
  out <- defaults
  out$m_vals <- validate_nonempty_vector(parse_integer_vector_arg(args$m_vals, defaults$m_vals), "m_vals")
  out$T_grid <- validate_nonempty_vector(parse_numeric_vector_arg(args$T_grid, defaults$T_grid), "T_grid")
  out$dgp_types <- validate_supported_values(
    validate_nonempty_vector(parse_character_vector_arg(args$dgp_types, defaults$dgp_types), "dgp_types"),
    streaming_curve_supported_dgp_types(),
    "dgp_types"
  )
  out$alpha_levels <- validate_nonempty_vector(parse_numeric_vector_arg(args$alpha_levels, defaults$alpha_levels), "alpha_levels")
  out$include_alt_detectors <- parse_bool_arg(args$include_alt_detectors, defaults$include_alt_detectors)
  out$overwrite_existing <- parse_bool_arg(args$overwrite_existing, defaults$overwrite_existing)
  out$nsim <- parse_scalar_integer_arg(args$nsim, defaults$nsim)
  out$q_cap <- parse_scalar_integer_arg(args$q_cap, defaults$q_cap)
  out$batch_size <- parse_scalar_integer_arg(args$batch_size, defaults$batch_size)
  out$basis_k <- parse_scalar_integer_arg(args$basis_k, defaults$basis_k)
  out$page_length_grid_size <- parse_scalar_integer_arg(args$page_length_grid_size, defaults$page_length_grid_size)
  out$weighted_length_grid_size <- parse_scalar_integer_arg(args$weighted_length_grid_size, defaults$weighted_length_grid_size)
  out$finite_eval_grid_size <- parse_scalar_integer_arg(args$finite_eval_grid_size, defaults$finite_eval_grid_size)
  out$exact_page_scan <- parse_bool_arg(args$exact_page_scan, defaults$exact_page_scan)
  out$exact_weighted_scan <- parse_bool_arg(args$exact_weighted_scan, defaults$exact_weighted_scan)
  out
}

resolve_power_run_settings <- function(args) {
  scenario <- as.character(args$scenario %||% "level_shift")
  defaults <- default_power_run_settings(scenario = scenario)
  out <- defaults
  out$scenario <- match.arg(scenario, streaming_curve_supported_power_scenarios())
  out$m_vals <- validate_nonempty_vector(parse_integer_vector_arg(args$m_vals, defaults$m_vals), "m_vals")
  out$T_grid <- validate_nonempty_vector(parse_numeric_vector_arg(args$T_grid, defaults$T_grid), "T_grid")
  out$dgp_types <- validate_supported_values(
    validate_nonempty_vector(parse_character_vector_arg(args$dgp_types, defaults$dgp_types), "dgp_types"),
    streaming_curve_supported_dgp_types(),
    "dgp_types"
  )
  out$s_star_vals <- validate_nonempty_vector(parse_integer_vector_arg(args$s_star_vals, defaults$s_star_vals), "s_star_vals")
  out$alpha_levels <- validate_nonempty_vector(parse_numeric_vector_arg(args$alpha_levels, defaults$alpha_levels), "alpha_levels")
  out$include_alt_detectors <- parse_bool_arg(args$include_alt_detectors, defaults$include_alt_detectors)
  out$overwrite_existing <- parse_bool_arg(args$overwrite_existing, defaults$overwrite_existing)
  out$nsim <- parse_scalar_integer_arg(args$nsim, defaults$nsim)
  out$q_cap <- parse_scalar_integer_arg(args$q_cap, defaults$q_cap)
  out$batch_size <- parse_scalar_integer_arg(args$batch_size, defaults$batch_size)
  out$basis_k <- parse_scalar_integer_arg(args$basis_k, defaults$basis_k)
  out$page_length_grid_size <- parse_scalar_integer_arg(args$page_length_grid_size, defaults$page_length_grid_size)
  out$weighted_length_grid_size <- parse_scalar_integer_arg(args$weighted_length_grid_size, defaults$weighted_length_grid_size)
  out$finite_eval_grid_size <- parse_scalar_integer_arg(args$finite_eval_grid_size, defaults$finite_eval_grid_size)
  out$exact_page_scan <- parse_bool_arg(args$exact_page_scan, defaults$exact_page_scan)
  out$exact_weighted_scan <- parse_bool_arg(args$exact_weighted_scan, defaults$exact_weighted_scan)

  global_delta <- parse_numeric_vector_arg(args$delta_vals, default = NULL)
  delta_map <- stats::setNames(vector("list", length(out$dgp_types)), out$dgp_types)
  for (nm in out$dgp_types) {
    dgp_key <- paste0("delta_vals_", nm)
    default_delta <- defaults$delta_map[[nm]]
    delta_map[[nm]] <- if (!is.null(global_delta)) {
      global_delta
    } else {
      parse_numeric_vector_arg(args[[dgp_key]], default = default_delta)
    }
  }
  out$delta_map <- delta_map
  out
}

build_power_param_grid <- function(settings) {
  rows <- list()
  idx <- 1L
  for (dgp_type in settings$dgp_types) {
    delta_vals <- settings$delta_map[[dgp_type]]
    for (m in settings$m_vals) {
      for (T_value in settings$T_grid) {
        for (s_star in settings$s_star_vals) {
          for (delta in delta_vals) {
            rows[[idx]] <- data.frame(
              m = as.integer(m),
              T_value = as.numeric(T_value),
              dgp_type = dgp_type,
              s_star = as.integer(s_star),
              delta = as.numeric(delta),
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}
