rGumbelFrechet <- function(n, theta, alpha) {
  cop <- gumbelCopula(param = theta, dim = n)
  U <- as.vector(rCopula(1, cop))
  X <- qfrechet(U, shape = alpha)
  return(X)
}