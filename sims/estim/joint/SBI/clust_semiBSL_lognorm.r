# ==============================================================================
# 0. Libraries
# ==============================================================================
library(BSL)
library(stabledist)
library(parallel)
library(coda)
library(here)

# ==============================================================================
# 1. Definitions (Functions)
# ==============================================================================

# Gumbel Generator (Frailty Representation)
rGumbV <- function(n, theta) {
  if (theta < 1) return(NULL)
  val_gamma <- (cos(pi / (2 * theta)))^theta
  V <- rstable(n = 1, alpha = 1 / theta, beta = 1,
               gamma = val_gamma, delta = 0, pm = 1)
  if (V <= 0) V <- .Machine$double.eps
  E <- rexp(n)
  U <- exp(-(E / V)^(1 / theta))
  return(U)
}

# Simulator
fn_sim_gumbel <- function(theta, n_obs) {
  mu <- theta[1]; sigma <- theta[2]; theta_cop <- theta[3]
  U <- rGumbV(n_obs, theta_cop)
  if (is.null(U) || any(is.nan(U))) return(NA)
  X <- qlnorm(U, meanlog = mu, sdlog = sigma)
  return(X)
}

# Summary Statistics
fn_sum_stats <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  s_max <- (max(x) - m) / md
  return(c(m, md, s_max))
}

# Log Prior
fn_log_prior <- function(theta) {
  mu <- theta[1]; sigma <- theta[2]; cop_theta <- theta[3]
  if (sigma <= 0 || cop_theta <= 1) return(-Inf)
  lp_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)
  lp_sigma <- dgamma(sigma, shape = 2, scale = 2, log = TRUE)
  lp_theta <- dgamma(cop_theta - 1, shape = 2, scale = 2, log = TRUE)
  return(lp_mu + lp_sigma + lp_theta)
}

# ==============================================================================
# 2. Data & Bounds
# ==============================================================================
set.seed(123)

N_OBS <- 1000 
true_params <- c(0, 1, 2) 
y_obs <- fn_sim_gumbel(true_params, n_obs = N_OBS)

theta_bounds <- cbind(
  c(-Inf, 0, 1),
  c( Inf, Inf, Inf)
)

# ==============================================================================
# 3. PILOT RUN
# ==============================================================================
cat("Step 1: Pilot Run...\n")

model_base <- newModel(
  fnSim = fn_sim_gumbel,
  fnSum = fn_sum_stats,
  fnLogPrior = fn_log_prior,
  theta0 = c(0.5, 1.5, 1.5), 
  simArgs = list(n_obs = N_OBS),
  thetaNames = c("mu", "sigma", "theta_cop")
)

pilot_res <- bsl(
  y = y_obs,
  n = 300, 
  M = 1000, 
  model = model_base,
  covRandWalk = diag(c(0.01, 0.01, 0.01)),
  method = "semiBSL",
  logitTransformBound = theta_bounds,
  verbose = FALSE
)

acc_samples <- getTheta(pilot_res)
tune_samples <- acc_samples[-(1:200), ]
tuned_cov <- cov(tune_samples)

if (any(is.na(tuned_cov)) || det(tuned_cov) <= 0) {
  tuned_cov <- diag(c(0.005, 0.005, 0.005))
} else {
  tuned_cov <- tuned_cov * (2.38^2 / 3)
}

# ==============================================================================
# 4. PARALLEL PRODUCTION RUN WITH LOGGING
# ==============================================================================

start_points <- list(
  c( 0.1, 1.0, 1.1),
  c( 0.2, 1.2, 1.8),
  c( 0.0, 0.9, 2.1),
  c(-0.1, 1.1, 1.5)
)

# ---- Create Log Directory ----
work_dir <- "/home/ferreroa/work_rainstsimu/extrCopula"
sbi_dir <- here(work_dir, "sims", "estim", "joint", "SBI")

log_dir <- here(sbi_dir, "logs")
dir.create(log_dir, showWarnings = FALSE)

# ---- Setup Parallel Cluster ----
n_cores <- min(length(start_points), detectCores() - 1)
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(BSL)
  library(stabledist)
})

clusterExport(cl, varlist = c(
  "y_obs", "tuned_cov", "N_OBS", "theta_bounds",      
  "rGumbV", "fn_sim_gumbel", 
  "fn_sum_stats", "fn_log_prior"
))

# ---- Worker with Full Logging ----
run_chain_worker <- function(start_val, chain_id, log_dir) {
  
  log_file <- file.path(log_dir, paste0("chain_", chain_id, ".log"))
  
  sink(log_file, append = TRUE, split = TRUE)
  sink(log_file, append = TRUE, type = "message")
  
  cat("\n=================================================\n")
  cat("CHAIN", chain_id, "STARTED at", as.character(Sys.time()), "\n")
  cat("Initial value:", paste(start_val, collapse = ", "), "\n")
  cat("=================================================\n\n")
  
  start_time <- Sys.time()
  
  res <- tryCatch({
    
    local_model <- newModel(
      fnSim = fn_sim_gumbel,
      fnSum = fn_sum_stats,
      fnLogPrior = fn_log_prior,
      theta0 = start_val,
      simArgs = list(n_obs = N_OBS),
      thetaNames = c("mu", "sigma", "theta_cop")
    )
    
    bsl(
      y = y_obs,
      n = 500,
      M = 5000,
      model = local_model,
      covRandWalk = tuned_cov,
      method = "semiBSL",
      logitTransformBound = theta_bounds,
      verbose = TRUE
    )
    
  }, error = function(e) {
    cat("\n!!! ERROR in Chain", chain_id, ":\n")
    cat(conditionMessage(e), "\n")
    return(NULL)
  })
  
  end_time <- Sys.time()
  
  cat("\n=================================================\n")
  cat("CHAIN", chain_id, "FINISHED at", as.character(end_time), "\n")
  cat("Total runtime:",
      as.numeric(difftime(end_time, start_time, units = "mins")),
      "minutes\n")
  cat("=================================================\n\n")
  
  sink()
  sink(type = "message")
  
  return(res)
}

# ---- Run Parallel Chains ----
cat("Running production chains...\n")

results_list <- parLapply(
  cl,
  seq_along(start_points),
  function(i) {
    run_chain_worker(
      start_val = start_points[[i]],
      chain_id  = i,
      log_dir   = log_dir
    )
  }
)

stopCluster(cl)

# ==============================================================================
# 5. Diagnostics
# ==============================================================================
mcmc_list <- mcmc.list(lapply(results_list, function(res) {
  samps <- getTheta(res)
  mcmc(samps[-(1:1000), ])
}))

plot(mcmc_list)
