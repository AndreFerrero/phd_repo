#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(copula)
library(evd)
library(future)
library(future.apply)

set.seed(123)

# ============================
# Paths
# ============================
work_dir <- "/home/ferreroa/work_rainstsimu/extrCopula"
dist_dir <- file.path(work_dir, "sims", "copula", "gumbel_distortion")
plots_dir <- file.path(dist_dir, "plots")
res_dir   <- file.path(dist_dir, "res")
common_dir <- file.path(work_dir, "sims", "common")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# Progress log file
log_file <- file.path(dist_dir, "progress.log")

write("======== NEW RUN STARTED ========", file = log_file)
write(sprintf("Start time: %s", Sys.time()), file = log_file, append = TRUE)

# ============================
# Helper functions
# ============================
source(file.path(common_dir, "handy_funs.r"))  # provides rCopFrechet etc.

# ============================
# TEST Simulation settings
# ============================
n <- 20000     
B <- 1000
R <- 500      
alpha <- 2
scenarios <- c("iid", "theta_1.5", "theta_2.5")

# ============================
# Parallel plan (over repetitions)
# ============================
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
plan(multisession, workers = n_cores)  # robust on HPC clusters

write(sprintf("Using %d cores", n_cores), file = log_file, append = TRUE)

# ============================
# Container for results
# ============================
all_results_by_scenario <- list()

# ============================
# Simulation loop
# ============================
for (scen in scenarios) {

  write(sprintf("\n=== Scenario %s started at %s ===", scen, Sys.time()),
        file = log_file, append = TRUE)

  cat("Starting scenario:", scen, "...\n")

  scenario_results <- future_lapply(seq_len(R), function(r) {

    # ===== progress every 50 reps =====
    if (r %% 50 == 0) {
      msg <- sprintf("%s | %s | reached repetition %d / %d",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     scen, r, R)
      write(msg, file = log_file, append = TRUE)
    }

    # ==================================
    # Generate maxima
    # ==================================
    max_vec <- numeric(B)

    if (scen == "iid") {
      for (i in seq_len(B)) {
        X <- rfrechet(n, shape = alpha)
        max_vec[i] <- max(X)
      }

    } else {
      tval <- as.numeric(sub("theta_", "", scen))
      for (i in seq_len(B)) {
        Gcop <- gumbelCopula(param = tval, dim = n)
        X <- rCopFrechet(alpha, Gcop)
        max_vec[i] <- max(X)
      }
    }

    # GEV fit
    fit <- tryCatch(
      fgev(max_vec, std.err = FALSE)$estimate,
      error = function(e) c(NA, NA, NA)
    )

    data.frame(
      repetition = r,
      location = fit[1],
      scale    = fit[2],
      shape    = fit[3]
    )

  }, future.seed = TRUE)

  all_results_by_scenario[[scen]] <- bind_rows(scenario_results)

  write(sprintf("=== Scenario %s finished at %s ===", scen, Sys.time()),
        file = log_file, append = TRUE)
}

# ============================
# Combine and save
# ============================
all_results_long <- bind_rows(lapply(names(all_results_by_scenario), function(scen) {
  df <- all_results_by_scenario[[scen]]
  df$scenario <- scen
  df
}))

save(all_results_long, file = file.path(res_dir, "clust_gumbel_wholemax_res.Rdata"))

write("\n======== RUN FINISHED ========", file = log_file, append = TRUE)
write(sprintf("End time: %s", Sys.time()), file = log_file, append = TRUE)

# Shutdown workers
plan(sequential)

# ============================
# Quick plot + save
# ============================
all_long <- all_results_long %>%
  pivot_longer(cols = c(location, scale, shape),
               names_to = "parameter",
               values_to = "estimate")

# Remove outliers
clean_long <- all_long
clean_long$keep <- TRUE

params <- unique(clean_long$parameter)
for (p in params) {
  idx <- clean_long$parameter == p
  est <- clean_long$estimate[idx]
  q_low  <- quantile(est, 0.005, na.rm = TRUE)
  q_high <- quantile(est, 0.995, na.rm = TRUE)
  clean_long$keep[idx] <- est >= q_low & est <= q_high
}

clean_long <- clean_long[clean_long$keep, ]
clean_long$keep <- NULL

my_cols <- c("iid" = "#66c2a5", "theta_1.5" = "#fc8d62", "theta_2.5" = "#8da0cb")

p <- ggplot(clean_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Distribution of GEV parameter estimates (cluster test)",
    subtitle = sprintf("Test run: n = %d, B = %d, R = %d", n, B, R),
    x = NULL, y = "Estimated parameter"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

save(p, file = file.path(plots_dir, "sim_plot.Rdata"))
