 
required_packages <- c("fdapace", "parallel", "fda", "zoo", "sandwich", "MASS", "fastmatrix")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(fdapace)
library(parallel)
library(fda)
library(zoo)
library(MASS)
library(fastmatrix)
 

################################################################################
# Zhu et al. (2025) — Adjusted-Range Based KS Statistic
# Includes critical values for gamma ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}
  
# ============================================================================== 
# Critical Value Lookup Table
# ============================================================================== 
get_critical_value_rsms <- function(T.chan, alpha = 0.05, gamma = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.1, 0.10)) stop("alpha must be 0.05 or 0.10.")
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15.")
  if (d < 1 || d > 8) stop("Only dimensions d = 1 to 8 are supported.")
  
  # Format alpha to fixed-point character
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  # Column map: alpha → T.chan index
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  # Validate mapping
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  # Retrieve index
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamma == 0, 0, 8)
  
  # Critical value matrix (RSMS)
  crit_mat <- matrix(c(
    2.0,1.5,2.7,2.0,3.2,2.6,3.8,2.9,    2.8,1.9,3.3,2.3,4.1,2.7,4.2,3.3,
    3.3,2.5,5.1,4.0,10.5,8.2,21.4,15.4, 5.1,4.0,7.9,6.2,20.1,15.0,46.5,33.3,
    4.1,3.3,7.0,5.5,14.4,11.7,27.0,21.0, 5.9,4.8,9.4,7.8,24.1,19.2,52.3,41.3,
    4.8,4.0,7.6,6.4,16.3,13.2,31.8,26.0, 6.8,5.3,11.3,9.2,24.1,20.1,69.2,51.9,
    5.4,4.5,9.2,7.3,18.0,14.7,34.2,28.3, 7.8,6.3,12.3,10.0,29.7,24.9,68.7,55.7,
    6.1,5.0,9.9,8.1,21.6,17.6,39.6,33.4, 9.0,7.3,14.1,11.7,33.6,27.8,73.4,59.6,
    7.2,6.0,10.6,8.8,21.2,18.5,39.7,35.2, 9.5,8.1,15.9,13.1,35.1,29.3,76.6,62.1,
    7.2,6.1,11.8,10.3,24.7,21.2,46.1,39.3, 9.8,8.8,16.2,13.8,39.4,32.4,84.6,69.4
  ), nrow = 8, byrow = TRUE)
  
  return(crit_mat[d, col_idx + offset])
}

# ============================================================================== 
# Adjusted-Range KS Test Statistic Function
# ============================================================================== 
 

rsms.statistic.fpca <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {

  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
    # Dimensions
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)  # Number of components (used as 'd')
  
  # Step 1: Demean using training sample
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Range-based variance normalization
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8  # Avoid division by zero
  V <- m * diag(1 / (range_vals^2))
  
  # Step 4: Compute scaled KS matrix
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- KS_mat * scaling
  KS_stat <- max(KS_scaled)
  
  # Step 5: Compare to critical value
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamma = gamm, d = K)
  reject <- KS_stat > crit_val
  
  return(list(statistic = KS_stat, critical_value = crit_val, reject = reject))
}
 
# ============================================================================== 
rsms.statistic.fpca.alt <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {
  input.vec <- as.matrix(input.vec)  # ensure matrix
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Handle vector case separately
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m ])
    
    cumsum_train <- cumsum(input_demeaned[1:m])
    cumsum_monitor <- cumsum(input_demeaned[(m+1):(m + m*T.chan)])
    
    range_val <- max(cumsum_train) - min(cumsum_train)
    if (range_val < 1e-8) range_val <- 1e-8
    V_inv <- m / (range_val^2)
    
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
      ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
    
    KS_scaled <- (cumsum_monitor^2) * V_inv * scaling
    KS_stat <- max(KS_scaled)
    
    crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamma = gamm, d = 1)
    reject <- KS_stat > crit_val
    first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
    
    return(list(
      statistic = KS_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # Multivariate case
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  cumsum_train <- apply(input_demeaned[1:m, , drop = FALSE], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE], 2, cumsum)
  
  cumsum_train <- t(cumsum_train)       # m x K
  cumsum_monitor <- t(cumsum_monitor)   # m*T.chan x K
  
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8
  V <- m * diag(1 / (range_vals^2))
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- diag(KS_mat) * scaling
  KS_stat <- max(KS_scaled)
  
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamma = gamm, d = K)
  reject <- KS_stat > crit_val
  first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
  
  return(list(
    statistic = KS_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

rsms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {
 
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
   # Dimensions
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training data
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Range-based variance normalization
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8
  V <- m * diag(1 / (range_vals^2))
  
  # Step 4: Compute test sequence over monitoring period
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  # Apply test
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- diag(KS_mat) * scaling
  KS_stat <- max(KS_scaled)
  
  # Step 5: Rejection rule
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamma = gamm, d = K)
  reject <- KS_stat > crit_val
  
  # Step 6: First time of rejection
  first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
  
  return(list(
    statistic = KS_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}


################################################################################
# Chan et al.'s (2020) Self-normalized statistics 
# Includes critical values for gamma ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}
 

# ============================================================================== 
# Critical Value Lookup Table
# ============================================================================== 

get_critical_value_ssms <- function(T.chan, alpha = 0.05, gamma = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.1, 0.10)) stop("alpha must be 0.05 or 0.10.")
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15.")
  if (d < 1 || d > 8) stop("Only dimensions d = 1 to 8 are supported.")
  
  # Format alpha for lookup
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamma == 0, 0, 8)
  
  crit_mat <- matrix(c(
    32.6,21.1,44.6,32.5,47.2,34.7,57.3,41.0,  50.2,32.9,51.2,37.0,64.6,44.7,74.2,50.5,
    78.0,50.6,95.1,67.3,123.9,89.4,128.2,95.9, 91.3,64.9,107.8,81.6,122.7,96.4,124.7,97.0,
    110.1,83.4,143.5,116.6,186.5,137.2,208.6,163.2, 138.8,106.7,185.4,137.6,223.3,161.6,222.0,171.5,
    151.3,120.9,204.5,165.8,263.4,208.4,303.6,247.8, 208.7,164.4,246.4,191.5,296.5,227.7,334.7,267.2,
    208.3,167.3,279.8,227.1,350.6,277.2,367.8,306.9, 272.2,209.6,341.8,263.0,382.5,313.1,398.2,328.7,
    290.4,221.4,340.8,291.9,451.2,377.8,493.0,406.8, 357.6,285.8,430.2,350.2,518.5,422.5,550.8,459.2,
    330.3,269.0,458.6,385.4,558.9,471.1,626.3,518.8, 459.3,356.6,543.7,444.3,611.7,522.2,636.8,525.7,
    383.7,325.7,572.2,486.2,686.5,582.6,751.7,607.8, 509.9,430.9,642.6,543.2,759.0,625.3,808.3,668.5
  ), nrow = 8, byrow = TRUE)
  
  return(crit_mat[d, col_idx + offset])
}

 
#get_critical_value_ssms(T.chan = 5, alpha = 0.05, gamma = 0.15, d = 4)
 

#get_critical_value_ssms(T.chan = 2, alpha = 0.10, gamma = 0, d = 6)
 
# ============================================================================== 
ssms.statistic.fpca  <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) { 
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sum over monitoring period (sn.k)
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)  # (m*T.chan) x K
  
  # Step 3: Compute denominator matrix D (based on training sample)
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  
  # Step 4: Compute the M-statistic vector
  library(MASS)
  Mt.mat <- sn.k %*% ginv(D) %*% t(sn.k)  # size (m*T.chan) x (m*T.chan)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2)
  Mt.vec <- diag(Mt.mat) * scaling
  
  # Step 5: Final statistic
  Mt_stat <- max(Mt.vec)
  
  # Step 6: Compare to critical value from SSMS table
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- Mt_stat > crit_val
  
  return(list(statistic = Mt_stat, critical_value = crit_val, reject = reject))
} 

# ==============================================================================  

ssms.statistic.fpca.alt <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) {
  input.vec <- as.matrix(input.vec)
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # ---- Univariate case ----
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m])
    
    # Step 2: cumulative sums of monitoring sample
    sn.k <- cumsum(input_demeaned[(m+1):(m + m*T.chan), ])
    
    # Step 3: Denominator matrix D
    cumsum_train <- cumsum(input_demeaned[1:m])
    D <- m^(-2) * sum(cumsum_train^2)
    
    # Step 4: Compute test sequence
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2)
    Mt.vec <- (sn.k^2 / D) * scaling
    
    # Step 5–7: Compute and return
    Mt_stat <- max(Mt.vec)
    crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = 1)
    reject <- Mt_stat > crit_val
    first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
    
    return(list(
      statistic = Mt_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # ---- Multivariate case ----
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE], 2, cumsum)
  cumsum_train <- apply(input_demeaned[1:m, , drop = FALSE], 2, cumsum)
  
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  Mt.mat <- sn.k %*% MASS::ginv(D) %*% t(sn.k)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2)
  Mt.vec <- diag(Mt.mat) * scaling
  
  Mt_stat <- max(Mt.vec)
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- Mt_stat > crit_val
  first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
  
  return(list(
    statistic = Mt_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}



ssms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) { 
 input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums of monitoring sample
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Denominator matrix D
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  
  # Step 4: Compute test sequence
  Mt.mat <- sn.k %*% MASS::ginv(D) %*% t(sn.k)
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2)
  Mt.vec <- diag(Mt.mat) * scaling
  
  # Step 5: Max test statistic
  Mt_stat <- max(Mt.vec)
  
  # Step 6: Critical value and rejection
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- Mt_stat > crit_val
  
  # Step 7: First rejection time
  first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
  
  return(list(
    statistic = Mt_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

# get_critical_value_rsms(T.chan = 5, alpha = 0.05, gamma = 0.15, d = 4)


# get_critical_value_rsms(T.chan = 2, alpha = 0.10, gamma = 0, d = 6)







################################################################################
# Standard HAC methods 
# Includes critical values for gamma ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}

# ==============================================================================
# Critical Value Lookup Table for CSMS  
# ==============================================================================

get_critical_value_csms <- function(T.chan, alpha = 0.05, gamma = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.10, 0.1)) stop("alpha must be 0.05 or 0.10.")
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15.")
  if (d < 1 || d > 8) stop("Only dimensions d = 1 to 8 are supported.")
  
  # Normalize alpha and T.chan keys
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  # Map alpha + T.chan → column index
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamma == 0, 0, 8)  # gamma = 0: cols 1–8; gamma = 0.15: cols 9–16
  
  # Critical value matrix for CSMS
  crit_mat <- matrix(c(
    2.2,1.7,3.4,2.7,3.8,3.0,4.9,3.6,   3.1,2.5,3.8,3.1,4.5,3.7,4.9,3.8,
    3.8,3.1,5.6,4.5,12.0,9.8,23.4,19.4, 5.3,4.3,9.0,7.1,18.0,14.9,32.8,26.9,
    4.4,3.8,7.5,6.1,15.9,13.5,28.9,24.2, 6.6,5.5,10.6,9.0,21.7,18.1,38.8,32.5,
    6.0,5.2,8.8,7.4,18.3,15.7,35.7,29.7, 7.7,6.4,12.3,10.8,26.6,21.9,47.6,40.8,
    6.7,5.6,10.2,8.6,21.0,18.0,37.5,32.6, 8.8,7.8,14.0,11.6,27.7,24.4,54.7,44.7,
    7.2,6.2,11.3,9.7,23.6,20.2,42.9,37.2, 10.0,8.5,15.3,13.2,31.0,26.8,62.4,53.0,
    7.9,7.1,12.7,11.1,26.8,23.2,50.5,43.9, 11.2,9.8,17.2,15.1,36.9,31.1,66.3,55.9,
    8.8,7.9,13.7,12.0,29.7,25.8,54.8,46.3, 11.9,10.2,18.7,16.2,38.6,32.7,70.8,63.8
  ), nrow = 8, byrow = TRUE)
  
  return(crit_mat[d, col_idx + offset])
}


# get_critical_value_csms(T.chan = 5, alpha = 0.05, gamma = 0.15, d = 4)


# get_critical_value_csms(T.chan = 2, alpha = 0.10, gamma = 0, d = 6)
# ==============================================================================

# ==============================================================================
# CSMS Test Statistic with HAC Normalization (based on training sample)
# ==============================================================================

csms.statistic.fpca <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) {
  library(sandwich)  # For HAC estimator
  library(MASS)      # For ginv()
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training sample mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Estimate HAC long-run variance (Newey-West via dummy regression)
  training_sample <- input_demeaned[1:m, ]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)  # Intercept-only model for each PC
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  # Step 3: Cumulative sum of monitoring period
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), ]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  
  # Step 4: Build CSMS statistic
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamma)
  
  # Final CSMS statistic
  csms_values <- diag(cumsum_monitor %*% ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  # Step 5: Compare with critical value from CSMS table
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- csms_stat > crit_val
  
  return(list(statistic = csms_stat, critical_value = crit_val, reject = reject))
}

# ==============================================================================
csms.statistic.fpca.alt <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) {
  library(sandwich)  # HAC estimator
  library(MASS)      # Generalized inverse
  
  input.vec <- as.matrix(input.vec)
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # === Univariate case ===
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m])
    
    # Step 2: Estimate HAC variance using Newey-West (variance only)
    ts_sample <- ts(input_demeaned[1:m, 1])
    lm_fit <- lm(ts_sample ~ 1)
    hac_var <- m * as.numeric(NeweyWest(lm_fit, prewhite = FALSE))
    
    if (is.na(hac_var) || is.nan(hac_var) || hac_var < 1e-8) {
      stop("HAC variance estimation failed or is too small.")
    }
    
    # Step 3: Cumulative sum of monitoring period
    cumsum_monitor <- cumsum(input_demeaned[(m+1):(m + m*T.chan), 1])
    
    # Step 4: Scaling
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
      ((time_idx / m) / (1 + time_idx / m))^(-2 * gamma)
    
    csms_values <- (cumsum_monitor^2 / hac_var) * scaling
    csms_stat <- max(csms_values)
    
    # Step 5: Critical value
    crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = 1)
    reject <- csms_stat > crit_val
    first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
    
    return(list(
      statistic = csms_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # === Multivariate case ===
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  training_sample <- input_demeaned[1:m, , drop = FALSE]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  cumsum_monitor <- t(cumsum_monitor)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamma)
  
  csms_values <- diag(cumsum_monitor %*% MASS::ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- csms_stat > crit_val
  first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
  
  return(list(
    statistic = csms_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

csms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, alpha = 0.05, gamma = 0) {
  library(sandwich)  # For HAC estimator
  library(MASS)      # For ginv()
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training sample mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Estimate HAC long-run variance (Newey-West via dummy regression)
  training_sample <- input_demeaned[1:m, ]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)  # Intercept-only model for each PC
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  # Step 3: Cumulative sum of monitoring period
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), ]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  
  # Step 4: Build CSMS statistic sequence
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamma)
  
  csms_values <- diag(cumsum_monitor %*% ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  # Step 5: Compare with critical value
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamma = gamma, d = K)
  reject <- csms_stat > crit_val
  first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
  
  return(list(
    statistic = csms_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

 

