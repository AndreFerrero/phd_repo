# =====================================================
# Simplified Gumbel Copula MLE using copula::fitCopula
# =====================================================

library(copula)

set.seed(123)

# ----------------------------
# 1. Simulation parameters
# ----------------------------
B <- 20             # number of blocks
n <- 200              # dimension of each block
theta_true <- 2     # true parameter

# ----------------------------
# 2. Simulate data
# ----------------------------
cop <- gumbelCopula(theta_true, dim = n)
U <- rCopula(B, cop)   # B i.i.d. samples of dimension n

# U is a B × n matrix, suitable for fitCopula()

# ----------------------------
# 3. Fit copula using built-in MLE
# ----------------------------
fit <- fitCopula(
  copula = gumbelCopula(dim = n),
  data   = U,
  method = "ml"
)

# ----------------------------
# 4. Output results
# ----------------------------
theta_hat <- coef(fit)

cat("Estimated Gumbel parameter:", theta_hat, "\n")
cat("True parameter:            ", theta_true, "\n")
