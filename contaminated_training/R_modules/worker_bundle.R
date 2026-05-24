# ==============================================================
# worker_bundle.R -- source all helper files on a PSOCK worker
# ============================================================== 

ROOT_W <- if (exists("ROOT2", inherits = TRUE)) ROOT2 else dirname(dirname(normalizePath(getwd(), winslash = "/", mustWork = FALSE)))
source(file.path(ROOT_W, "R", "project_paths.R"), local = FALSE)
source(file.path(ROOT_W, "R", "dependencies.R"), local = FALSE)
source(file.path(ROOT_W, "R", "utils.R"), local = FALSE)
source(file.path(ROOT_W, "R", "critical_values_lookup.R"), local = FALSE)
source(file.path(ROOT_W, "R", "fpca_pipeline.R"), local = FALSE)
source(file.path(ROOT_W, "R", "method_catalog.R"), local = FALSE)
source(file.path(ROOT_W, "R", "scenarios.R"), local = FALSE)
source(file.path(ROOT_W, "R", "contamination.R"), local = FALSE)
source(file.path(ROOT_W, "R", "detectors_main.R"), local = FALSE)
source(file.path(ROOT_W, "R", "detectors_alt.R"), local = FALSE)
source(file.path(ROOT_W, "R", "simulation_core.R"), local = FALSE)
source(file.path(ROOT_W, "R", "genData.R"), local = FALSE)
NULL
