# ==============================================================
# preflight_checks.R -- setup validation before long simulation runs
# ==============================================================

preflight_worker_bootstrap <- function(root, nworkers = 2L) {
  nworkers <- as.integer(max(1L, nworkers))
  if (nworkers <= 1L) return(invisible(TRUE))

  tryCatch({
    cl <- make_psock_cluster(nworkers)
    on.exit(stop_psock_cluster(cl), add = TRUE)
    ROOT2 <- root
    parallel::clusterExport(cl, varlist = c("ROOT2"), envir = environment())
    parallel::clusterEvalQ(cl, {
      source(file.path(ROOT2, "R", "worker_bundle.R"), local = FALSE)
      TRUE
    })
    TRUE
  }, error = function(e) {
    warning("PSOCK worker bootstrap failed; simulation will fall back to serial execution if needed: ",
            conditionMessage(e), call. = FALSE)
    FALSE
  })
}

run_simulation_preflight <- function(root,
                                     settings,
                                     main_catalog,
                                     alt_catalog = data.frame(),
                                     ncores = 1L,
                                     check_worker_bootstrap = TRUE) {
  ensure_simulation_packages(install_if_missing = FALSE)

  cv_main <- load_main_critical_values(root)
  validate_main_critical_value_coverage(
    cv_main = cv_main,
    main_catalog = main_catalog,
    T_grid = settings$T_grid,
    q_grid = seq_len(as.integer(settings$q_cap)),
    alpha_levels = settings$alpha_levels
  )

  cv_alt <- NULL
  alt_warning <- NULL
  if (!is.null(alt_catalog) && nrow(alt_catalog) > 0L) {
    alt_result <- tryCatch({
      cv_alt_local <- load_alt_critical_values(root)
      validate_alt_critical_value_coverage(
        cv_alt = cv_alt_local,
        alt_catalog = alt_catalog,
        T_grid = settings$T_grid,
        q_grid = seq_len(as.integer(settings$q_cap)),
        alpha_levels = settings$alpha_levels
      )
      list(cv_alt = cv_alt_local, warning = NULL)
    }, error = function(e) {
      msg <- paste0(
        "Benchmark detectors are being disabled for this run because their critical-value file is unavailable or invalid. ",
        "Regenerate it with `Rscript generate_alt_detector_critical_values.R`. Details: ",
        conditionMessage(e)
      )
      warning(msg, call. = FALSE)
      list(cv_alt = NULL, warning = msg)
    })
    cv_alt <- alt_result$cv_alt
    alt_warning <- alt_result$warning
  }

  worker_bootstrap_ok <- TRUE
  if (isTRUE(check_worker_bootstrap) && as.integer(ncores) > 1L) {
    worker_bootstrap_ok <- isTRUE(preflight_worker_bootstrap(root, nworkers = min(2L, as.integer(ncores))))
  }

  list(
    cv_main = cv_main,
    cv_alt = cv_alt,
    alt_warning = alt_warning,
    worker_bootstrap_ok = worker_bootstrap_ok
  )
}
