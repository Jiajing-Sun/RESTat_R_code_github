# ==============================================================
# alt_detector_weights.R -- detector-specific weights
# ============================================================== 

normalize_omega_name <- function(weight) {
  w <- toupper(trimws(as.character(weight)[1L]))
  if (w %in% c("U", "UNIFORM", "NONE", "CONST", "CONSTANT", "ONE", "1")) return("Uniform")
  if (w %in% c("INVSQRT", "INV_SQRT", "SQRT")) return("InvSqrt")
  if (w %in% c("INVLIN", "INV_LINEAR", "INVLINEAR", "LIN")) return("InvLinear")
  stop("Unknown weighted-CUSUM omega. Use one of: Uniform, InvSqrt, InvLinear.")
}

make_weighted_cusum_omega <- function(length_s, weight = "Uniform", eps = 1e-8) {
  w <- normalize_omega_name(weight)
  L <- pmax(as.numeric(length_s), eps)
  if (w == "Uniform") return(rep(1, length(L)))
  if (w == "InvSqrt") return(1 / sqrt(L))
  if (w == "InvLinear") return(1 / L)
  stop("Unknown weighted-CUSUM omega.")
}

normalize_scale_weight_name <- function(weight) {
  w <- toupper(trimws(as.character(weight)[1L]))
  if (w %in% c("EQUAL", "CONST", "CONSTANT", "UNIFORM", "ONE", "1")) return("Equal")
  if (w %in% c("INVSQRTH", "INV_SQRT_H", "SQRTH")) return("InvSqrtH")
  stop("Unknown multiscale weight. Use one of: Equal, InvSqrtH.")
}

make_multiscale_weight <- function(h, weight = "Equal", eps = 1e-8) {
  w <- normalize_scale_weight_name(weight)
  h <- pmax(as.numeric(h), eps)
  if (w == "Equal") return(rep(1, length(h)))
  if (w == "InvSqrtH") return(1 / sqrt(h))
  stop("Unknown multiscale weight.")
}
