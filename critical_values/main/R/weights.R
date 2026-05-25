# ==============================================================
# weights.R -- CvM weights for finite and open-end monitoring
# ============================================================== 

normalize_weight_name <- function(weight) {
  w <- toupper(trimws(as.character(weight)[1L]))
  if (w %in% c("U", "UNIFORM", "CONST", "CONSTANT", "ONE", "WU", "W_U", "1")) return("U")
  if (w %in% c("EARLY", "W_EARLY")) return("Early")
  if (w %in% c("MID", "W_MID")) return("Mid")
  if (w %in% c("LATE", "W_LATE")) return("Late")
  stop("Unknown weight. Use one of: U, Early, Mid, Late.")
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

make_cvm_weight_open_x <- function(x, weight = "U") {
  w <- normalize_weight_name(weight)

  if (w == "U") return(rep(1, length(x)))
  if (w == "Early") return(2 * (1 - x))
  if (w == "Late") return(2 * x)
  if (w == "Mid") return(6 * x * (1 - x))

  stop("Unknown open-end weight.")
}

make_cvm_weight_open_s <- function(s, weight = "U") {
  w <- normalize_weight_name(weight)

  if (w == "U") return((1 + s)^(-2))
  if (w == "Early") return(2 * (1 + s)^(-3))
  if (w == "Late") return(2 * s * (1 + s)^(-3))
  if (w == "Mid") return(6 * s * (1 + s)^(-4))

  stop("Unknown open-end weight.")
}
