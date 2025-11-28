rCopFrechet <- function(alpha, cop) {
  # Given a copula, return Frechet(alpha) sample
  U <- as.vector(rCopula(1, cop))
  X <- qfrechet(U, shape = alpha)
  return(X)
}

rGumbV <- function(n, theta) {
    require(stabledist)

    # Gumbel V r.v.
    V <- rstable(
        n = 1,
        alpha = 1/theta,
        beta = 1,
        gamma = cospi(1/(2 * theta))^theta,
        delta = 0,
        pm = 1
    )

    E <- rexp(n, V) 

    # Gumbel generator
    gum_psi <- function(t, theta){
        exp(- t ^ (1/theta))
    }

    #(U_1, ..., U_n) ~ C_psi
    U <- gum_psi(E, theta)

    return(U)
}

block_max <- function(x, block_size) {
  # return block maxima given the block size
  sapply(split(x, ceiling(seq_along(x) / block_size)), max)
}