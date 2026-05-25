# ==============================================================
# paper_baseline.R -- structured manuscript baseline for audit
# ============================================================== 

streaming_curve_paper_dirname <- function() {
  "20260328_SUN_ZHU_Haer_LIN_Online_Monitoring_via_Streaming_Curves (13)"
}

paper_baseline_spec <- function(root) {
  paper_dir <- file.path(root, streaming_curve_paper_dirname())
  list(
    sources = list(
      manuscript_tex = file.path(paper_dir, "Online_Monitoring_via_Streaming_Curves.tex"),
      benchmark_tex = file.path(paper_dir, "benchmark_detector_critical_values_streaming.tex"),
      appendix_tex = file.path(paper_dir, "appendix_critical_values_streaming.tex")
    ),
    null_design = list(
      dgp_types = c("BB", "fIID", "fMA1"),
      m_vals = c(100L, 200L, 500L, 1000L),
      T_grid = c(1, 2, 5, 10),
      gamma_vec = c(0, 0.15),
      alpha_nominal = 0.05,
      nsim = 1000L
    ),
    power_design = list(
      dgp_types = c("BB", "fIID", "fMA1"),
      m_vals = 500L,
      T_grid = c(1, 2),
      s_star_vals = c(50L, 200L),
      scenarios = c("level_shift", "smooth_change", "abrupt_local_change", "sinusoidal_change"),
      basis_dimension = streaming_curve_basis_dimension(),
      localized_basis_index = localized_change_basis_index_default(),
      nsim = 1000L,
      delta_map = list(
        level_shift = list(
          BB = c(0.050, 0.058, 0.065, 0.073, 0.080, 0.088, 0.095, 0.103, 0.110, 0.118, 0.125),
          fIID = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029),
          fMA1 = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029)
        ),
        smooth_change = list(
          BB = c(0.050, 0.058, 0.065, 0.073, 0.080, 0.088, 0.095, 0.103, 0.110, 0.118, 0.125),
          fIID = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029),
          fMA1 = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029)
        ),
        abrupt_local_change = list(
          BB = c(0.100, 0.400, 0.700, 1.000, 1.300, 1.600),
          fIID = c(0.05, 0.08, 0.11, 0.14, 0.17, 0.20),
          fMA1 = c(0.05, 0.08, 0.11, 0.14, 0.17, 0.20)
        ),
        sinusoidal_change = list(
          BB = c(0.050, 0.058, 0.065, 0.073, 0.080, 0.088, 0.095, 0.103, 0.110, 0.118, 0.125),
          fIID = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029),
          fMA1 = c(0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029)
        )
      ),
      formulas = list(
        smooth_change = "Delta * (t - t*) / (mT)",
        sinusoidal_change = "Delta * sin(pi * (t - t*) / (mT))"
      )
    ),
    critical_values = list(
      q_grid = 1:30,
      T_grid = c(1, 2, 5, 10),
      alpha_levels = c(0.01, 0.05, 0.10),
      gamma_vec = c(0, 0.15),
      finite_nrep = 10000L,
      open_nrep = 5000L,
      finite_kmax = 10000L,
      open_train_grid_size = 1500L,
      open_monitor_grid_size = 2000L,
      benchmark_weighted_omega = "InvSqrt",
      mosum_h_vec = c(0.10, 0.20),
      multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
      multiscale_scale_names = c("Equal")
    ),
    manuscript_notes = list(
      hac_bandwidth_text = "Silverman's Rule of Thumb",
      expected_arl_artifact_prefix = "ARL_output_"
    )
  )
}
