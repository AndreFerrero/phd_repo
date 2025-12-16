source("libs/packages.R")

build_logposterior <- function(copula, margin, param_map, data,
                               transform = NULL,
                               inverse_transform = NULL,
                               log_jacobian = NULL) {
  
  function(param_in) {
    
    # Apply inverse transform if given
    param <- if (!is.null(inverse_transform)) inverse_transform(param_in) else param_in
    
    param_m <- param[param_map$margin]
    param_c <- param[param_map$copula]
    
    # 1. Prior
    lprior <- copula$log_prior(param_c) +
          margin$log_prior(param_m)
    
    if (!is.finite(lprior)) return(-Inf)
    
    # 2. Likelihood (copula + marginal)
    u <- margin$cdf(data, param_m)
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)  # numerical stability
    
    llik <- sum(margin$log_density(data, param_m)) +
          copula$log_density(u, param_c)
    
    # 3. Jacobian adjustment if in transformed space
    jac <- if (!is.null(log_jacobian)) log_jacobian(param_in) else 0
    
    lprior + llik + jac
  }
}
