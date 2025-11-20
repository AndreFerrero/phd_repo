library(copula)

set.seed(123)
B <- 20       # number of blocks
n <- 5        # block size
theta_true <- 2

# Simulate block-dependent data
simulate_blocks <- function(B, n, theta) {
  cop <- gumbelCopula(theta, dim=n)
  lapply(1:B, function(b) as.numeric(rCopula(1, cop)))
}

blocks <- simulate_blocks(B, n, theta_true)

# ----------------------------
# Closed-form diagonal MLE
# ----------------------------
Y <- sapply(blocks, max)  # maxima of each block
theta_hat_closed <- log(n) / (log(B) - log(sum(-log(Y))))

cat("True theta:", theta_true, "\n")
cat("Closed-form DMLE:", theta_hat_closed, "\n")
