source("libs/packages.R")


build_simulator <- function(copula, margin, param_map) {

  function(param, n) {

    param_m <- param[param_map$margin]
    param_c <- param[param_map$copula]

    U <- copula$simulate_u(param_c, n)
    if (is.null(U)) return(NA)

    X <- margin$simulate(U, param_m)
    if (any(!is.finite(X))) return(NA)

    X
  }
}
