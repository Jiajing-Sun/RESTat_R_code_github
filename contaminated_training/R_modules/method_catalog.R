# ==============================================================
# method_catalog.R -- metadata for all simulation methods
# ============================================================== 

build_main_method_catalog <- function(gamma_vec = c(0, 0.15),
                                      cvm_weights = c("U", "Early", "Mid", "Late"),
                                      standardizers = c("HAC", "SSMS", "RSMS")) {
  rows <- list(); idx <- 1L
  for (std in standardizers) {
    for (g in gamma_vec) {
      rows[[idx]] <- data.frame(family = "Main", standardizer = std, detector = "Main", type = "KS",
                                gamma = as.numeric(g), weight_name = "", bandwidth_h = NA_real_, omega_name = "",
                                hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
    for (w in cvm_weights) {
      rows[[idx]] <- data.frame(family = "Main", standardizer = std, detector = "Main", type = "CvM",
                                gamma = NA_real_, weight_name = normalize_weight_name(w), bandwidth_h = NA_real_, omega_name = "",
                                hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$method_label <- mapply(method_label_main, out$standardizer, out$type, out$gamma, out$weight_name, USE.NAMES = FALSE)
  out$method_group <- mapply(method_group_from_row, out$family, out$type, out$detector, USE.NAMES = FALSE)
  out$method_id <- sanitize_tag(tolower(gsub("[^A-Za-z0-9]+", "_", out$method_label)))
  out
}

build_alt_method_catalog <- function(gamma_vec = c(0, 0.15),
                                     standardizers = c("HAC", "SSMS", "RSMS"),
                                     mosum_h_vec = c(0.10, 0.20),
                                     weighted_omega_names = c("InvSqrt"),
                                     multiscale_h_sets = list(H050_100_200 = c(0.05, 0.10, 0.20)),
                                     multiscale_scale_names = c("Equal")) {
  rows <- list(); idx <- 1L
  for (std in standardizers) {
    for (g in gamma_vec) {
      rows[[idx]] <- data.frame(family = "Benchmark", standardizer = std, detector = "FullCUSUM", type = "KS",
                                gamma = as.numeric(g), weight_name = "", bandwidth_h = NA_real_, omega_name = "",
                                hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(family = "Benchmark", standardizer = std, detector = "PageCUSUM", type = "KS",
                                gamma = as.numeric(g), weight_name = "", bandwidth_h = NA_real_, omega_name = "",
                                hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
    for (h in mosum_h_vec) {
      rows[[idx]] <- data.frame(family = "Benchmark", standardizer = std, detector = "MOSUM", type = "KS",
                                gamma = NA_real_, weight_name = "", bandwidth_h = as.numeric(h), omega_name = "",
                                hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
    for (g in gamma_vec) {
      for (om in weighted_omega_names) {
        rows[[idx]] <- data.frame(family = "Benchmark", standardizer = std, detector = "WeightedCUSUM", type = "KS",
                                  gamma = as.numeric(g), weight_name = "", bandwidth_h = NA_real_, omega_name = normalize_omega_name(om),
                                  hset_name = "", scale_weight_name = "", stringsAsFactors = FALSE); idx <- idx + 1L
      }
    }
    for (hs in names(multiscale_h_sets)) {
      for (sw in multiscale_scale_names) {
        rows[[idx]] <- data.frame(family = "Benchmark", standardizer = std, detector = "MultiscaleMOSUM", type = "KS",
                                  gamma = NA_real_, weight_name = "", bandwidth_h = NA_real_, omega_name = "",
                                  hset_name = hs, scale_weight_name = normalize_scale_weight_name(sw), stringsAsFactors = FALSE); idx <- idx + 1L
      }
    }
  }
  out <- do.call(rbind, rows)
  out$method_label <- mapply(method_label_alt, out$standardizer, out$detector, out$gamma, out$bandwidth_h,
                             out$omega_name, out$hset_name, out$scale_weight_name, USE.NAMES = FALSE)
  out$method_group <- mapply(method_group_from_row, out$family, out$type, out$detector, USE.NAMES = FALSE)
  out$method_id <- sanitize_tag(tolower(gsub("[^A-Za-z0-9]+", "_", out$method_label)))
  out
}
