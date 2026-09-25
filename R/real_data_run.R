# =============================================================================
# real_data_run.R
# Fits the GVEM-based IRT model to the real dataset (ARC benchmark).
#
# NOTE: The indices i and j are used in the reverse order compared to the
# paper. Specifically, i indexes items in the paper but indexes groups here,
# and j indexes groups in the paper but indexes items here.
#
# USAGE:
#   1. Place arc_resp.RData and dummy_matrix.RData in the data/ folder.
#   2. Source main.R before running this script, or ensure it is sourced below.
#   3. Outputs are saved to results/real_data/.
# =============================================================================
 
rm(list = ls())
options(scipen = 999)

library(mirt)
library(dplyr)
library(tidyr)

source("R/main.R")

# -----------------------------------------------------------------------------
# Load real data
# -----------------------------------------------------------------------------
resp         <- readRDS("data/arc_resp.RData")
dummy_matrix <- readRDS("data/dummy_matrix.RData")

# Dimensions for reference
I <- 5
J <- 100
S <- 6005

dim(dummy_matrix)

# -----------------------------------------------------------------------------
# Initialization
# -----------------------------------------------------------------------------
group     <- as.factor(resp$groupfull)
all_start <- start_value(resp, J = 100)

# -----------------------------------------------------------------------------
# Fit GVEM (single lambda)
# -----------------------------------------------------------------------------
GVEM_out <- GVEM_main(
  resp          = resp,
  all_start     = all_start,
  dummy_matrix  = dummy_matrix,
  lambda        = 30,
  c             = 0.04 * S,
  iter_criteria = 5e2,
  tau_criteria  = 1e-3,
  rho_N         = 1e-2,
  rho_N2        = 0.1,
  debias        = FALSE
)

# -----------------------------------------------------------------------------
# Inspect results
# -----------------------------------------------------------------------------
GVEM_out$sig2_b_j
GVEM_out$sig2_theta_i
GVEM_out$theta_i

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
main_alpha <- as.data.frame(t(GVEM_out$alpha_main))

write.csv(main_alpha, "results/real_data/main_alpha_results.csv")