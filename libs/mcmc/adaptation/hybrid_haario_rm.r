#' Hybrid adaptive scheme: Robbins–Monro + Haario covariance
#'
#' @param target_accept Target acceptance probability for Robbins–Monro
#' @param gamma0 Initial Robbins–Monro step size
#' @param kappa Robbins–Monro decay exponent (> 0.5)
#' @param eps Small jitter added to covariance diagonal
#' @param t0 Start adapting empirical covariance (Haario)
#'
#' @return Adaptation object with update(state, param, accept, iter)
adapt_hybrid <- function(
  target_accept = 0.234,
  gamma0 = 1,
  kappa = 1,
  eps = 1e-6,
  t0 = 1000
) {

  list(
    update = function(state, param, accept, iter) {

      d <- length(param)

      # --- Robbins–Monro: adapt scale ---
      if (iter > 1) {
        gamma_t <- gamma0 / iter^kappa
        log_scale <- log(state$scale) + gamma_t * (as.numeric(accept) - target_accept)
        state$scale <- exp(log_scale)
      }

      # --- Haario: empirical covariance adaptation ---
      if (iter > t0) {
        delta <- param - state$mean
        state$mean <- state$mean + delta / state$t
        state$cov  <- state$cov + (tcrossprod(delta) - state$cov) / state$t
        state$t <- state$t + 1
      }

      # --- Update proposal covariance ---
      # Apply Robbins–Monro scale to empirical covariance (with jitter)
      state$Sigma <- state$scale^2 * (state$cov + eps * diag(d))

      state
    }
  )
}
