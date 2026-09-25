# =============================================================================
# simulation.R
# Simulation study for TriEval.
#
# NOTE: The indices i and j are used in the reverse order compared to the
# paper. Specifically, i indexes items in the paper but indexes groups here,
# and j indexes groups in the paper but indexes items here.
#
# USAGE:
#   1. Place main_alpha_results.csv in the data/ folder.
#   2. Outputs are written to results/simulation/.
#   3. Run this script from the repository root.
# =============================================================================
 
rm(list = ls())

library(mirt)
library(dplyr)
library(tidyr)
library(parallel)

source("R/main.R")

# -----------------------------------------------------------------------------
# Load true alpha parameters
# -----------------------------------------------------------------------------
main_alpha      <- read.csv("data/main_alpha_results.csv")
main_alpha_real <- as.matrix(main_alpha[, 2:6])


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Expand and Perturb a Parameter Matrix
#'
#' Repeats columns of a matrix (cycling 1, 2, ..., 5, 1, 2, ...) to reach
#' \code{target_cols} columns, then adds small Gaussian noise to all entries.
#'
#' @param mat        Numeric matrix. Source parameter matrix.
#' @param rows        Integer. Number of rows to use from \code{mat}.
#' @param target_cols Integer. Desired number of columns in the output.
#' @param noise_sd    Numeric. Standard deviation of the added noise. Default 0.1.
#'
#' @return A numeric matrix of dimension \code{rows x target_cols}.
expand_cols <- function(mat, rows, target_cols, noise_sd = 0.1) {
  col_indices <- rep_len(1:5, target_cols)
  new_mat     <- mat[1:rows, col_indices]
  n_rows      <- nrow(new_mat)
  noise       <- matrix(rnorm(n_rows * target_cols, mean = 0, sd = noise_sd),
                        nrow = n_rows,
                        ncol = target_cols)
  return(new_mat + noise)
}


#' Simulate Prompt-Effect IRT Data
#'
#' Generates a synthetic dataset under the model with
#' model- and task-level random effects driven by a factorial prompt design.
#'
#' @param levels   Integer vector of length 5. Number of levels for each
#'                 prompt factor (V1--V5).
#' @param I        Integer. Number of groups (models).
#' @param sig_beta Numeric. Standard deviation of the item intercept prior.
#' @param rr       Integer. Random seed.
#'
#' @return A named list with:
#'   \item{resp}{Data frame of item responses with a \code{groupfull} column.}
#'   \item{dummy_matrix}{Design matrix of prompt factor dummies (S x P).}
#'   \item{S}{Integer. Total number of prompt conditions.}
#'   \item{sim_true}{Named list of true parameter values used in simulation.}
sim_data_prompt <- function(levels, I, sig_beta, rr) {
  
  set.seed(rr)
  
  # Build full factorial design matrix
  df <- expand.grid(
    V1 = factor(1:levels[1]),
    V2 = factor(1:levels[2]),
    V3 = factor(1:levels[3]),
    V4 = factor(1:levels[4]),
    V5 = factor(1:levels[5])
  )
  
  dummy_matrix       <- model.matrix(~ V1 + V2 + V3 + V4 + V5, data = df)[, -1]
  dummy_matrix_tilde <- cbind(1, dummy_matrix)
  
  S <- nrow(dummy_matrix)
  J <- 100
  P <- ncol(dummy_matrix)
  
  # ---------------------------------------------------------------------------
  # True parameter generation
  # ---------------------------------------------------------------------------
  sim_theta_i <- c(rep(c(-1, 1), I / 2))
  sim_sigma_i <- c(rep(c(0, 1), each = I / 2))
  
  sim_alpha       <- matrix(NA, nrow = I, ncol = (P + 1))
  sim_alpha[, 1]  <- sim_theta_i
  sim_alpha[, -1] <- t(expand_cols(mat       = main_alpha_real,
                                   rows      = P,
                                   target_cols = I,
                                   noise_sd  = 0.1))
  
  sim_mu_theta <- sim_alpha %*% t(dummy_matrix_tilde)
  
  # Vectorized generation of prompt-specific model abilities
  sim_theta_is          <- t(sim_mu_theta +
                               sim_sigma_i * matrix(rnorm(I * S), nrow = I, ncol = S))
  colnames(sim_theta_is) <- paste0("Model_", 1:I)
  
  sim_d      <- rnorm(J, 0, sig_beta)
  sim_sigma_j <- c(rep(0, J / 2), rep(1, J / 2))
  sim_b_js   <- t(sim_d + sim_sigma_j * matrix(rnorm(J * S), nrow = J, ncol = S))
  
  # ---------------------------------------------------------------------------
  # Response generation (vectorized)
  # ---------------------------------------------------------------------------
  theta_long <- as.vector(sim_theta_is)
  b_long     <- sim_b_js[rep(1:S, times = I), ]
  logits     <- theta_long - b_long
  probs      <- 1 / (1 + exp(-logits))
  resp_raw   <- matrix(as.integer(runif(length(probs)) < probs),
                       nrow = S * I, ncol = J)
  
  # ---------------------------------------------------------------------------
  # Format output
  # ---------------------------------------------------------------------------
  resp             <- as.data.frame(resp_raw)
  colnames(resp)   <- paste0("Item_", 1:J)
  rownames(resp)   <- paste0("Model_", rep(1:I, each = S),
                             "_s",    rep(1:S, times = I))
  resp$groupfull   <- rep(1:I, each = S)
  
  sim_true <- list(
    sim_theta_i = sim_theta_i,
    sim_sigma_i = sim_sigma_i,
    sim_alpha   = sim_alpha,
    sim_sigma_j = sim_sigma_j,
    sig_beta    = sig_beta
  )
  
  return(list(
    resp         = resp,
    dummy_matrix = dummy_matrix,
    S            = S,
    sim_true     = sim_true
  ))
}


# -----------------------------------------------------------------------------
# Evaluation Metrics
# -----------------------------------------------------------------------------

#' Compute Estimation Error Metrics
#'
#' @param x    Numeric vector or matrix of estimated values.
#' @param true Numeric scalar, vector, or matrix of true values.
#'
#' @return A named list with RMSE, bias, and mean absolute bias.
eva_metrics <- function(x, true) {
  rmse  <- sqrt(mean(as.vector((x - true)^2)))
  bias  <- mean(as.vector((x - true)))
  abias <- mean(as.vector(abs(x - true)))
  return(list(rmse = rmse, bias = bias, abs = abias))
}


# -----------------------------------------------------------------------------
# Simulation Configuration
# -----------------------------------------------------------------------------
configs <- expand.grid(
  levels_idx = 1:2,
  I          = c(6, 20),
  sig_beta   = c(1, sqrt(2))
)

level_options <- list(
  c(10, 5, 5, 5, 2),
  c(5,  5, 2, 2, 2)
)

# Lambda grids per configuration
all_lambda_list <- list(
  seq(20, 30, 1),
  seq(1,  5,  0.5),
  seq(20, 60, 4),
  seq(1,  10, 1),
  seq(20, 30, 1),
  seq(1,  5,  0.5),
  seq(20, 60, 4),
  seq(1,  10, 1)
)

# Output directory — update this path before running
main_dir <- "results/simulation"


# -----------------------------------------------------------------------------
# Main Simulation Loop
# -----------------------------------------------------------------------------
for (n in 1:nrow(configs)) {
  for (rep in 1:50) {
    
    conf     <- configs[n, ]
    levs     <- level_options[[conf$levels_idx]]
    S_val    <- prod(levs)
    
    folder_name <- sprintf("I%d_sig%.1f_S%d", conf$I, conf$sig_beta, S_val)
    folder_path <- file.path(main_dir, folder_name)
    
    if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    
    setwd(folder_path)
    
    cat(sprintf("Running: I=%d, sig_beta=%.1f, S=%d, rep=%d\n",
                conf$I, conf$sig_beta, prod(levs), rep))
    
    # Simulate data
    sim_data <- sim_data_prompt(levels   = levs,
                                I        = conf$I,
                                sig_beta = conf$sig_beta,
                                rr       = 1000 * n + rep)
    
    saveRDS(sim_data$sim_true, paste0("sim_true_rep", rep, ".RData"))
    
    group     <- as.factor(sim_data$resp$groupfull)
    all_start <- start_value(sim_data$resp, J = 100)
    
    ncore     <- 10
    S         <- sim_data$S
    all_lambda <- all_lambda_list[[n]]
    
    # Fit RGVEM across lambda grid in parallel
    RGVEM_returns <- mclapply(all_lambda, function(lambda_val) {
      RGVEM_main(
        resp           = sim_data$resp,
        all_start      = all_start,
        dummy_matrix   = sim_data$dummy_matrix,
        lambda         = lambda_val,
        c              = 0.1 * S,
        iter_criteria  = 5e2,
        tau_criteria   = 1e-3,
        rho_N          = 1e-2,
        rho_N2         = 0.1,
        debias         = FALSE
      )
    }, mc.cores = ncore)
    
    # Select best model by GIC and BIC
    GIC1         <- sapply(RGVEM_returns, function(x) x$GIC)
    BIC1         <- sapply(RGVEM_returns, function(x) x$BIC)
    RGVEM_minGIC  <- RGVEM_returns[[which.min(GIC1)]]
    RGVEM_minBIC  <- RGVEM_returns[[which.min(BIC1)]]
    
    saveRDS(RGVEM_returns, paste0("RGVEM_returns_rep", rep, ".RData"))
    saveRDS(RGVEM_minGIC,  paste0("minGIC_rep",       rep, ".RData"))
    saveRDS(RGVEM_minBIC,  paste0("minBIC_rep",       rep, ".RData"))
  }
}


# -----------------------------------------------------------------------------
# Results Summary Loop
# -----------------------------------------------------------------------------
for (n in 1:nrow(configs)) {
  conf        <- configs[n, ]
  levs        <- level_options[[conf$levels_idx]]
  S_val       <- prod(levs)
  
  folder_name <- sprintf("I%d_sig%.1f_S%d", conf$I, conf$sig_beta, S_val)
  folder_path <- file.path(main_dir, folder_name)
  
  setwd(folder_path)
  
  sim_true <- readRDS(paste0("sim_true_rep", rep, ".RData"))
  RGVEM     <- readRDS(paste0("minGIC_rep",   rep, ".RData"))
  
  RGVEM$sig2_b_j
  RGVEM$sig2_theta_i
  
  eva_metrics(RGVEM$theta_i,    sim_theta_i)
  eva_metrics(RGVEM$alpha_main, 0.2)
}