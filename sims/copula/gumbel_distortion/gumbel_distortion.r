library(parallel)

set.seed(123)

plots_dir <- here("sims", "copula", "gumbel_distortion", "plots")
res_dir <- here("sims", "copula", "gumbel_distortion", "res")
common_dir <- here("sims", "common")

# ============================================================
# Helper functions
# ============================================================

source(here(common_dir, "handy_funs.r"))

# ============================================================
# Simulation settings
# ============================================================

n <- 20000
B <- 1000
block_size <- 200
(n_blocks <- n / block_size) 
alpha <- 2
scenarios <- c("iid", "theta_1.5", "theta_2.5")

# ============================================================
# Parallel setup
# ============================================================

cl <- makeCluster(parallel::detectCores() - 1)
clusterEvalQ(cl, {
  library(copula)
  library(evd)
})
clusterExport(cl, c("n", "alpha", "block_size", "block_max", "rCopFrechet"))

# ============================================================
# Simulation loop
# ============================================================

results_list <- list()

for (scen in scenarios) {
  cat("Running scenario:", scen, "...\n")

  if (scen == "iid") {
    res <- parSapply(cl, 1:B, function(i, n, alpha, block_size) {
      X <- rfrechet(n, shape = alpha)
      M <- block_max(X, block_size)
      fit <- tryCatch(fgev(M, std.err = FALSE)$estimate,
        error = function(e) c(NA, NA, NA)
      )
      fit
    }, n, alpha, block_size)
  } else {
    tval <- as.numeric(sub("theta_", "", scen))
    clusterExport(cl, c("tval"), envir = environment())
    clusterEvalQ(cl, {
      Gcop <- gumbelCopula(param = tval, dim = n)
    })

    res <- parSapply(cl, 1:B, function(i, alpha, block_size) {
      X <- rCopFrechet(alpha, Gcop)
      M <- block_max(X, block_size)
      fit <- tryCatch(fgev(M, std.err = FALSE)$estimate,
        error = function(e) c(NA, NA, NA)
      )
      fit
    }, alpha, block_size)
  }

  results_list[[scen]] <- t(res)
}

stopCluster(cl)

# ============================================================
# Combine results into tidy data frame
# ============================================================

all_results <- bind_rows(lapply(names(results_list), function(scen) {
  df <- as.data.frame(results_list[[scen]])
  colnames(df) <- c("location", "scale", "shape")
  df$scenario <- scen
  df
}))

# save(all_results, file = here(res_dir, "gumbel_res.Rdata"))
# load(file = here(res_dir, "gumbel_res.Rdata"))
# ============================================================
# Reshape and plot all parameters
# ============================================================

# Convert to long format
all_long <- all_results %>%
  pivot_longer(
    cols = c(location, scale, shape),
    names_to = "parameter",
    values_to = "estimate"
  )

# Color palette
my_cols <- c(
  "iid" = "#66c2a5",
  "theta_1.5" = "#fc8d62",
  "theta_2.5" = "#8da0cb"
)

# Plot faceted boxplots
ggplot(all_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.7, outlier.color = "grey40") +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Comparison of GEV Parameter Estimates - F ~ Frechet(2)",
    subtitle = paste("IID vs Gumbel", "-", "n =", n, "B =", B, "Block size =", block_size, "Number of blocks =", n_blocks),
    x = NULL, y = "Estimated parameter value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold")
  )

ggplot(all_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_violin(alpha = 0.7) +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Comparison of GEV Parameter Estimates - F ~ Frechet(2)",
    subtitle = paste(
      "IID vs Gumbel --",
      "Simulations =", B, " Observations =", n, " Block size =", block_size,
      " Number of blocks =", n / block_size
    ),
    x = NULL, y = "Estimated parameter value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold")
  )

