margin_lognormal <- list(

  name = "lognormal",

  simulate = function(u, param) {
    qlnorm(u, param["mu"], param["sigma"])
  },

  log_density = function(x, param) {
    dlnorm(x, param["mu"], param["sigma"], log = TRUE)
  },

  cdf = function(x, param) {
    plnorm(x, param["mu"], param["sigma"])
  },

  log_prior = function(param) {
    if (param["sigma"] <= 0) return(-Inf)
    dnorm(param["mu"], 0, 10, log = TRUE) +
      dlnorm(param["sigma"], 0, 1, log = TRUE)
  }
)
