
# Load packages
source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")

# Load builders
source("libs/models/builders/simulator.R")
source("libs/models/builders/logposterior.R")
source("libs/models/builders/logposterior_indep.R")

# Load MCMC machinery
source("libs/mcmc/run_chain.R")
source("libs/mcmc/run_parallel_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")
source("libs/mcmc/adaptation/robbins_monro.R")

# Simulate data
param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

set.seed(123)
n_obs <- 150

true_param_indep <- c(mu = 0, sigma = 1)
X_indep <- rlnorm(n_obs, true_param_indep[1], true_param_indep[2])

true_param_dep <- c(mu = 0, sigma = 1, theta = 2)
X_dep <- simulator(true_param_dep, n_obs)

# =============================================================================
# 2. Build log-posterior transformations
# =============================================================================

# --- FOR DEPENDENT MODEL (3 params) ---
g_dep <- function(param) {
  c(param["mu"], log(param["sigma"]), log(param["theta"] - 1))
}

g_inv_dep <- function(phi) {
  param <- c(mu = phi[1], sigma = exp(phi[2]), theta = exp(phi[3]) + 1)
  names(param) <- c("mu", "sigma", "theta")
  param
}

log_jacobian_dep <- function(phi) phi[2] + phi[3]

# --- FOR INDEPENDENT MODEL (2 params) ---
g_indep <- function(param) {
  c(param["mu"], log(param["sigma"])) # Only mu and sigma
}

g_inv_indep <- function(phi) {
  param <- c(mu = phi[1], sigma = exp(phi[2]))
  names(param) <- c("mu", "sigma")
  param
}

log_jacobian_indep <- function(phi) phi[2]

logpost_indep <- build_logposterior_indep(
  margin = margin_lognormal,
  param_map = param_map,
  data = X_indep,
  inverse_transform = g_inv_indep,
  log_jacobian = log_jacobian_indep
)

logpost_dep <- build_logposterior(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map,
  data = X_dep,
  inverse_transform = g_inv_dep,
  log_jacobian = log_jacobian_dep
)

phi_indep <- g_indep(true_param_indep)
phi_dep <- g_dep(true_param_dep)

proposal_indep <- proposal_gaussian_rw(Sigma = diag(0.01, 2))
proposal_dep <- proposal_gaussian_rw(Sigma = diag(0.1, 3))

n_mcmc <- 15000

# Number of chains
n_chains <- 4

# Overdispersed initial points
inits_indep <- lapply(1:n_chains, function(i) {
  phi_indep + rnorm(length(phi_indep), 0, 0.2)
})

inits_dep <- lapply(1:n_chains, function(i) {
  phi_dep + rnorm(length(phi_dep), 0, 0.2)
})

# Run parallel MCMC
res_par_indep <- run_parallel_chains(
  log_target  = logpost_indep,
  init_values = inits_indep,
  n_iter      = 10000,
  proposal    = proposal_indep,
  adapt       = adapt_haario(),
  transform   = g_inv_indep,
  n_cores     = n_chains,
  export      = c(
    "g_inv_indep",
    "copula_gumbel",
    "margin_lognormal",
    "X_indep",
    "log_jacobian_indep",
    "mh_step"
  )
)

res_par_dep <- run_parallel_chains(
  log_target  = logpost_dep,
  init_values = inits_dep,
  n_iter      = 10000,
  proposal    = proposal_dep,
  adapt       = adapt_haario(),
  transform   = g_inv_dep,
  n_cores     = n_chains,
  export      = c(
    "g_inv_dep",
    "copula_gumbel",
    "margin_lognormal",
    "X_dep",
    "log_jacobian_dep",
    "mh_step"
  )
)

# res_indep <- run_chain(
#   log_target = logpost_indep,
#   init = phi_indep,
#   n_iter = n_mcmc,
#   proposal = proposal_indep,
#   adapt = adapt_haario()
# )

# res_dep <- run_chain(
#   log_target = logpost_dep,
#   init = phi_dep,
#   n_iter = n_mcmc,
#   proposal = proposal_dep,
#   adapt = adapt_haario()
# )

full_joint_lik <- here("sims", "estim", "joint", "full_lik_bayes")
res_full_joint_lik <- here(full_joint_lik, "res")

# save(res_dep, res_indep, file = here(res_full_joint_lik, "indep_vs_dep_10kiter_haario_lognorm_gumbel.Rdata"))

samples_indep <- t(apply(res_indep$samples, 1, g_inv_indep))
samples_dep <- t(apply(res_dep$samples, 1, g_inv_dep))

mcmc_samples_indep <- mcmc(samples_indep)
mcmc_samples_dep <- mcmc(samples_dep)

burn_in <- nrow(samples_indep)/2
# burn_in <- 0
mcmc_clean_indep <- window(mcmc_samples_indep, start = burn_in + 1, thin = 1)
mcmc_clean_dep <- window(mcmc_samples_dep, start = burn_in + 1, thin = 1)

# Arrange 1 row and 3 columns
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  # smaller margins for tighter plots

plot(as.matrix(mcmc_clean_dep[, "mu"]), type = "l",
     ylab = "mu", xlab = "Iteration", main = "Trace of mu")

plot(as.matrix(mcmc_clean_dep[, "sigma"]), type = "l",
     ylab = "sigma", xlab = "Iteration", main = "Trace of sigma")

plot(as.matrix(mcmc_clean_dep[, "theta"]), type = "l",
     ylab = "theta", xlab = "Iteration", main = "Trace of theta")

# Reset par
par(mfrow = c(1,1))

# Arrange 1 row and 3 columns
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))  # smaller margins for tighter plots

plot(as.matrix(mcmc_clean_indep[, "mu"]), type = "l",
     ylab = "mu", xlab = "Iteration", main = "Trace of mu")

plot(as.matrix(mcmc_clean_indep[, "sigma"]), type = "l",
     ylab = "sigma", xlab = "Iteration", main = "Trace of sigma")

# Reset par
par(mfrow = c(1,1))


summary(mcmc_samples_indep)

# =============================================================================
# 4. Recover Limit Distribution G
# =============================================================================

# Define grid and helper function
y_grid <- seq(0.01, 100, length.out = 200)

compute_G <- function(params, grid, n) {
  # params: vector with mu, sigma, and optionally theta
  mu    <- params["mu"]
  sigma <- params["sigma"]
  theta <- if("theta" %in% names(params)) params["theta"] else 1
  
  F_y   <- plnorm(grid, meanlog = mu, sdlog = sigma)
  F_y   <- pmin(pmax(F_y, 1e-15), 1 - 1e-15)
  
  # Use your copula's diagonal function
  # Note: your copula_gumbel$diag expects (u, theta) and uses length(u) as n
  # We override it here to use the correct sample size n_obs
  F_y^(n^(1/theta))
}

# --- Calculate G for all posterior draws ---
G_dep   <- t(apply(samples_dep, 1, compute_G, grid = y_grid, n = n_obs))
G_indep <- t(apply(samples_indep, 1, compute_G, grid = y_grid, n = n_obs))

# --- True G for comparison ---
G_true_dep <- compute_G(true_param_dep, y_grid, n_obs)
G_true_indep <- compute_G(true_param_indep, y_grid, n_obs)

# =============================================================================
# 5. Visualization
# =============================================================================

summarize_G <- function(G_mat) {
  list(
    mean   = colMeans(G_mat),
    lower  = apply(G_mat, 2, quantile, 0.05),
    upper  = apply(G_mat, 2, quantile, 0.95)
  )
}

stats_dep   <- summarize_G(G_dep)
stats_indep <- summarize_G(G_indep)

par(mfrow = c(1, 2))
# Theta = 2
plot(y_grid, G_true_dep, type = "l", lwd = 2, lty = 2, col = "red",
     ylim = c(0, 1), xlab = "y", ylab = expression(P(M[n] <= y)),
     main = "Maxima Distribution G, theta = 2")

# Plot Dependent Model Posterior
lines(y_grid, stats_dep$mean, col = "blue", lwd = 2)
polygon(c(y_grid, rev(y_grid)), c(stats_dep$lower, rev(stats_dep$upper)),
        col = rgb(0, 0, 1, 0.1), border = NA)

# Theta = 1
plot(y_grid, G_true_indep, type = "l", lwd = 2, lty = 2, col = "red",
     ylim = c(0, 1), xlab = "y", ylab = expression(P(M[n] <= y)),
     main = "Maxima Distribution G, theta = 1")

# Plot Independent Model Posterior
lines(y_grid, stats_indep$mean, col = "darkgreen", lwd = 2)

polygon(c(y_grid, rev(y_grid)), c(stats_indep$lower, rev(stats_indep$upper)),
        col = rgb(0, 0, 1, 0.1), border = NA)
par(mfrow = c(1, 1))


