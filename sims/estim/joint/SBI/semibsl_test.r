# =============================================================================
# Test script: small MCMC run with Gumbel copula and Lognormal margin
# =============================================================================

# Load packages
source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")

# Load builders
source("libs/models/builders/simulator.R")
source("libs/models/builders/bsl_logposterior.R")
source("libs/models/builders/synthetic_semibsl.r")
source("libs/models/builders/synthetic_gaussian_bsl.r")

# Load MCMC machinery
source("libs/mcmc/run_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")
source("libs/mcmc/adaptation/robbins_monro.R")
source("libs/mcmc/adaptation/hybrid_haario_rm.R")

# =============================================================================
# 1. Simulate data using the simulator
# =============================================================================
param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

set.seed(123)
true_param <- c(mu = 0, sigma = 1, theta = 2)
n_obs <- 1000
X <- simulator(true_param, n_obs)

# =============================================================================
# 2. Build log-posterior
# =============================================================================
g <- function(param) {
  c(param["mu"], log(param["sigma"]), log(param["theta"] - 1))
}

g_inv <- function(phi) {
  param <- c(
    mu    = phi[1],
    sigma = exp(phi[2]),
    theta = exp(phi[3]) + 1
  )
  # Ensure names are exactly what we want
  names(param) <- c("mu", "sigma", "theta")
  param
}

med_mad_max <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  smax <- (max(x) - m) / md
  raw_max <- max(x)
  c(m, md, raw_max)
}

log_jacobian <- function(phi) phi[2] + phi[3]

logpost <- build_bsl_logposterior(
  copula = copula_gumbel, margin = margin_lognormal,
  param_map = param_map,
  data = X,
  simulator = simulator,
  sum_stats = med_mad_max,
  n_sim = 150,
  synthetic_loglik = synthetic_semibsl,
  inverse_transform = g_inv,
  log_jacobian = log_jacobian
)

# =============================================================================
# 3. Define proposal and run chain
# =============================================================================
phi_init <- g(c(mu = 0, sigma = 1, theta = 2))
p <- 3
Sigma0 <- matrix(0.1, p, p)
diag(Sigma0) <- 1
proposal <- proposal_gaussian_rw(Sigma0 = Sigma0)

res <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 10000,
  proposal = proposal,
  adapt = adapt_haario()
)

sbi_dir <- here("sims", "estim", "joint", "SBI")
sbi_res_dir <- here(sbi_dir, "res")

# save(res, file = here(sbi_res_dir, "semibsl_10kruns_150sims_1kobs_adapthaario.Rdata"))
load(here(sbi_res_dir, "semibsl_10kruns_150sims_1kobs_adapthaario.Rdata"))

# Convert samples back to natural space
samples_natural <- t(apply(res$samples, 1, g_inv))

mcmc_samples <- mcmc(samples_natural)

burn_in <- nrow(samples_natural)/2
# burn_in <- 0
mcmc_clean <- window(mcmc_samples, start = burn_in + 1, thin = 1)

# =============================================================================
# 4. Quick summaries
# =============================================================================
res$conv

cat("Acceptance rate:", res$accept_rate, "\n")
cat("Posterior means:\n")
print(colMeans(as.matrix(mcmc_clean)))

effectiveSize(mcmc_clean)
# Traceplot

# Arrange 1 row and 3 columns
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  # smaller margins for tighter plots

plot(as.matrix(mcmc_clean[, "mu"]), type = "l",
     ylab = "mu", xlab = "Iteration", main = "Trace of mu")

plot(as.matrix(mcmc_clean[, "sigma"]), type = "l",
     ylab = "sigma", xlab = "Iteration", main = "Trace of sigma")

plot(as.matrix(mcmc_clean[, "theta"]), type = "l",
     ylab = "theta", xlab = "Iteration", main = "Trace of theta")

# Reset par
par(mfrow = c(1,1))

pairs(as.matrix(mcmc_clean))

acf(as.matrix(mcmc_clean))

# Densities
par(mfrow = c(1, 3))

mu_dens_est <- density(as.matrix(mcmc_clean)[, "mu"])
plot(mu_dens_est,
    main = "Posterior Density of Mu",
    lwd = 2, col = "blue", xlab = expression(mu)
)
abline(v = true_param[1], col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

sigma_dens_est <- density(as.matrix(mcmc_clean)[, "sigma"])
plot(sigma_dens_est,
    main = "Posterior Density of Sigma",
    lwd = 2, col = "blue", xlab = expression(sigma)
)
abline(v = true_param[2], col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

theta_dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(theta_dens_est,
    main = "Posterior Density of Theta",
    lwd = 2, col = "blue", xlab = expression(theta)
)
abline(v = true_param[3], col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

par(mfrow = c(1, 1))
