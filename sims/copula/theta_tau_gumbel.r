plots_folder <- here("sims", "copula", "plots")

set.seed(123)

# Sample size
n <- 500

# Define different theta values (corresponding to tau)
theta_values <- c(2, 3, 4, 6)
tau_values <- 1 - 1/theta_values

# Create an empty list to store data
datasets <- list()

# Generate datasets
for(i in 1:length(theta_values)){
  gumbel_cop <- gumbelCopula(param = theta_values[i], dim = 2)
  data <- rCopula(n, gumbel_cop)
  datasets[[i]] <- data
}

png(filename = here(plots_folder, "theta_tau_gumbel.png"))
# Plotting datasets
par(mfrow = c(2,2))  # 2x2 grid

for(i in 1:length(datasets)){
  plot(datasets[[i]], main = paste0("Gumbel Copula\nθ=", round(theta_values[i],2),
                                     ", τ=", round(tau_values[i],2)),
       xlab = "U1", ylab = "U2", pch = 19, col = rgb(0,0,1,0.5))
}

dev.off()
