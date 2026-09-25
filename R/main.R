# =============================================================================
# main.R
# Core functions for regularized Gaussian variational expectation–maximization 
# (RGVEM) estimation
# =============================================================================

# NOTE: The indices i and j are used in the reverse order compared to the
# paper. Specifically, i indexes task in the paper but indexes model here,
# and j indexes model in the paper but indexes task here.

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Compute Starting Values
#'
#' Initializes starting values for the RGVEM algorithm.
#'
#' @param resp  A data frame of task responses (rows = respondents, 
#'              columns = tasks indicators, 
#'              with the last column as model indicator).
#' @param J     Integer. Number of tasks.
#'
#' @return A named list with:
#'   \item{theta_i}{Initial baseline model capability (vector of length I).}
#'   \item{sig2_theta_i}{Initial prompt-by-model variances (vector of length I).}
#'   \item{sig2_b_j}{Initial prompt-by-task variances (vector of length J).}
start_value <- function(resp, J) {
  group <- as.factor(resp$groupfull)
  
  start_theta_i   <- rep(1,   n_distinct(group))
  start_sig2_i    <- rep(0.5, n_distinct(group))
  start_sig2_j    <- rep(1,   J)
  
  start_list <- list(
    theta_i      = start_theta_i,
    sig2_theta_i = start_sig2_i,
    sig2_b_j     = start_sig2_j
  )
  
  return(start_list)
}


#' Jaakkola-Jordan Local Variational Parameter Function
#'
#' Computes the eta function used in the variational lower bound,
#' as defined in Jaakkola & Jordan (2000).
#'
#' @param y A numeric vector or matrix of squared variational parameters.
#'
#' @return A numeric vector or matrix of the same dimension as \code{y}.
eta <- function(y) {
  x <- sqrt(y)
  return((1 / (2 * x)) * ((1 / (1 + exp(-x))) - 1 / 2))
}


# -----------------------------------------------------------------------------
# Main RGVEM Function
# -----------------------------------------------------------------------------

#' regularized Gaussian variational expectation–maximization (RGVEM) Estimation
#'
#' Fits a TriEval model using a regularized Gaussian variational expectation–maximization algorithm with log-penalty.
#'
#' @param resp           A data frame of task responses with the last column as model indicator.
#' @param all_start      A named list of starting values from \code{start_value()}.
#' @param lambda         Numeric. log-penalty parameter.
#' @param c              Numeric. Scaling constant for the GIC criterion.
#' @param dummy_matrix   A numeric matrix of group-level covariates (I x P).
#' @param iter_criteria  Integer. Maximum number of EM iterations. Default 500.
#' @param tau_criteria   Numeric. Convergence threshold. Default 1e-3.
#' @param rho_N          Numeric. Threshold for hard-zeroing small variances
#'                       after convergence (sparsity detection).
#' @param rho_N2         Numeric. Minimum variance floor for likelihood
#'                       computation.
#' @param debias         Logical. If FALSE (default), apply sparsity thresholding
#'                       after convergence.
#'
#' @return A named list containing:
#'   \item{sig2_beta}{Estimated variance of overall task intercept prior.}
#'   \item{sig2_theta_i}{Estimated model-level prompt-by-model variances (length I).}
#'   \item{sig2_b_j}{Estimated task-level prompt-by-task variances (length J).}
#'   \item{alpha}{Augmented coefficient matrix (I x (P+1)), including intercepts.}
#'   \item{alpha_main}{Main prompt effect on capability, only (I x P).}
#'   \item{theta_i}{Estimated model capability intercepts (length I).}
#'   \item{likelihood}{Variational lower bound at convergence.}
#'   \item{BIC}{Bayesian Information Criterion.}
#'   \item{GIC}{Generalized Information Criterion.}
#'   \item{iter}{Number of iterations until convergence.}
#'   \item{tau}{Final convergence criterion value.}
RGVEM_main <- function(resp, all_start,
                      lambda, c,
                      dummy_matrix,
                      iter_criteria = 5e2, tau_criteria = 1e-3,
                      rho_N, rho_N2,
                      debias = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Data preparation
  # ---------------------------------------------------------------------------
  Ysij <- split(resp %>% select(-starts_with("group")), f = resp$groupfull)
  Ysij <- lapply(Ysij, as.matrix)
  Ysj  <- Reduce("+", Ysij)
  
  Ns <- as.vector(unlist(lapply(Ysij, nrow)))
  I  <- length(unique(resp$groupfull))
  J  <- ncol(resp %>% select(-starts_with("group")))
  P  <- ncol(dummy_matrix)
  
  dummy_matrix_tilde <- cbind(1, dummy_matrix)
  
  # ---------------------------------------------------------------------------
  # Initialize variational parameters
  # ---------------------------------------------------------------------------
  ksi_sij2        <- list()
  sig2_beta       <- list()
  sig2_theta_i    <- list()
  sig2_b_j        <- list()
  alpha           <- list()
  
  q_sig2_beta_0j  <- list()
  q_mu_beta_0j    <- list()
  q_sig2_theta_is <- list()
  q_mu_theta_is   <- list()
  q_sig2_b_js     <- list()
  q_mu_b_js       <- list()
  
  ksi_sij2[[1]]     <- lapply(Ns, function(s) matrix(0.1, nrow = s, ncol = J))
  sig2_beta[[1]]    <- 0.5
  sig2_theta_i[[1]] <- all_start$sig2_theta_i
  sig2_b_j[[1]]     <- all_start$sig2_b_j
  
  alpha[[1]]        <- matrix(0, nrow = I, ncol = (P + 1))
  alpha[[1]][, 1]   <- all_start$theta_i
  alpha[[1]][1:I, -1] <- 0
  
  q_sig2_beta_0j[[1]] <- rep(1, J)
  q_mu_beta_0j[[1]]   <- rep(0, J)
  
  q_sig2_theta_is[[1]] <- list()
  q_mu_theta_is[[1]]   <- list()
  Ysi_list              <- list()
  
  for (i in 1:I) {
    q_sig2_theta_is[[1]][[i]] <- rep(1, Ns[i])
    q_mu_theta_is[[1]][[i]]   <- rep(0, Ns[i])
    Ysi_list[[i]]             <- rowSums(Ysij[[i]])
  }
  
  Ysi <- do.call(rbind, Ysi_list)
  
  q_sig2_b_js[[1]] <- matrix(1, nrow = J, ncol = S)
  q_mu_b_js[[1]]   <- matrix(0, nrow = J, ncol = S)
  
  # Sparsity index sets
  out_sig2_b0    <- which(sig2_b_j[[1]] != 0)
  in_sig2_b0     <- which(sig2_b_j[[1]] == 0)
  out_sig2_theta <- which(sig2_theta_i[[1]] != 0)
  in_sig2_theta  <- which(sig2_theta_i[[1]] == 0)
  
  iter <- 1
  tau  <- 1
  
  # ---------------------------------------------------------------------------
  # EM Iterations
  # ---------------------------------------------------------------------------
  while ((tau >= tau_criteria) && (iter <= iter_criteria)) {
    
    # -- Initialize E-step storage -------------------------------------------
    q_sig2_beta_0j[[iter + 1]]  <- rep(NA, J)
    q_mu_beta_0j[[iter + 1]]    <- rep(NA, J)
    q_sig2_theta_is[[iter + 1]] <- list()
    q_mu_theta_is[[iter + 1]]   <- list()
    q_sig2_b_js[[iter + 1]]     <- matrix(NA, nrow = J, ncol = S)
    q_mu_b_js[[iter + 1]]       <- matrix(NA, nrow = J, ncol = S)
    
    # -- E-step: compute variational expectations ----------------------------
    eta_sij2_temp <- list()
    eta_si2_temp  <- list()
    ksi_theta_temp <- list()
    ksi_b_temp    <- list()
    ksi_b_s_temp  <- list()
    
    for (i in 1:I) {
      eta_sij2_temp[[i]]  <- eta(ksi_sij2[[iter]][[i]])
      eta_si2_temp[[i]]   <- rowSums(eta_sij2_temp[[i]])
      ksi_theta_temp[[i]] <- q_mu_theta_is[[iter]][[i]] * eta_sij2_temp[[i]]
      ksi_b_temp[[i]]     <- q_mu_b_js[[iter]] * t(eta_sij2_temp[[i]])
      ksi_b_s_temp[[i]]   <- colSums(ksi_b_temp[[i]])
    }
    
    eta_sj2_temp <- Reduce("+", eta_sij2_temp)
    
    # Update q(beta_0j)
    temp1               <- S * sig2_beta[[iter]] + sig2_b_j[[iter]]
    q_sig2_beta_0j[[iter + 1]] <- (sig2_b_j[[iter]] * sig2_beta[[iter]]) / temp1
    q_mu_beta_0j[[iter + 1]]   <- rowSums(q_mu_b_js[[iter]]) * (sig2_beta[[iter]] / temp1)
    
    # Update q(b_js)
    sig2_b_j_temp <- rep(sig2_b_j[[iter]], each = nrow(eta_sj2_temp))
    temp2         <- 1 + 2 * eta_sj2_temp * sig2_b_j_temp
    
    q_sig2_b_js[[iter + 1]][out_sig2_b0, ] <- t(sig2_b_j_temp / temp2)[out_sig2_b0, ]
    q_sig2_b_js[[iter + 1]][in_sig2_b0, ]  <- 0
    
    q_mu_beta_0j_temp <- rep(q_mu_beta_0j[[iter]], each = nrow(eta_sj2_temp))
    temp3 <- 2 * Reduce("+", ksi_theta_temp)
    temp4 <- Ysj - I / 2 - temp3
    temp5 <- q_mu_beta_0j_temp - sig2_b_j_temp * temp4
    
    q_mu_b_js[[iter + 1]][out_sig2_b0, ] <- t(temp5 / temp2)[out_sig2_b0, ]
    q_mu_b_js[[iter + 1]][in_sig2_b0, ]  <- q_mu_beta_0j[[iter]][in_sig2_b0]
    
    # Update q(theta_is)
    alpha_x_temp <- alpha[[iter]] %*% t(dummy_matrix_tilde)
    temp6        <- 1 + 2 * sig2_theta_i[[iter]] * do.call(rbind, eta_si2_temp)
    temp7        <- do.call(rbind, ksi_b_s_temp)
    temp8        <- Ysi - J / 2 + 2 * temp7
    temp9        <- alpha_x_temp + sig2_theta_i[[iter]] * temp8
    
    q_sig2_theta_is_temp <- sig2_theta_i[[iter]] / temp6
    q_mu_theta_is_temp   <- temp9 / temp6
    
    q_sig2_theta_is_temp[in_sig2_theta, ] <- 0
    q_mu_theta_is_temp[in_sig2_theta, ]   <- alpha_x_temp[in_sig2_theta, ]
    
    for (i in 1:I) {
      q_sig2_theta_is[[iter + 1]][[i]] <- q_sig2_theta_is_temp[i, ]
      q_mu_theta_is[[iter + 1]][[i]]   <- q_mu_theta_is_temp[i, ]
    }
    
    # -- M-step: update hyperparameters --------------------------------------
    ksi_sij2[[iter + 1]]     <- list()
    alpha[[iter + 1]]        <- alpha[[iter]]
    sig2_theta_i[[iter + 1]] <- sig2_theta_i[[iter]]
    sig2_b_j[[iter + 1]]     <- sig2_b_j[[iter]]
    sig2_beta[[iter + 1]]    <- sig2_beta[[iter]]
    
    tau_M <- 1
    
    while (tau_M > tau_criteria) {
      
      alpha_old        <- alpha[[iter + 1]]
      sig2_theta_i_old <- sig2_theta_i[[iter + 1]]
      sig2_b_j_old     <- sig2_b_j[[iter + 1]]
      sig2_beta_old    <- sig2_beta[[iter + 1]]
      
      m_temp1 <- q_mu_b_js[[iter + 1]]^2 + q_sig2_b_js[[iter + 1]]
      m_temp2 <- list()
      m_temp3 <- list()
      m_temp4 <- list()
      
      for (i in 1:I) {
        m_temp2[[i]]          <- q_mu_theta_is[[iter + 1]][[i]]^2 + q_sig2_theta_is[[iter + 1]][[i]]
        m_temp3[[i]]          <- sweep(q_mu_b_js[[iter + 1]], 2, q_mu_theta_is[[iter + 1]][[i]], "*")
        m_temp4[[i]]          <- sweep(-2 * m_temp3[[i]], 2, m_temp2[[i]], "+")
        ksi_sij2[[iter + 1]][[i]] <- t(m_temp4[[i]] + m_temp1)
      }
      
      # Update alpha (both the baseline model capability and model-level prompt effect)
      q_mu_theta_is_temp    <- do.call(rbind, q_mu_theta_is[[iter + 1]])
      alpha[[iter + 1]]     <- (q_mu_theta_is_temp %*% dummy_matrix_tilde) %*%
        solve(crossprod(dummy_matrix_tilde))
      
      # Update sig2_theta_i (with log penalty)
      m_temp5              <- (q_mu_theta_is_temp - alpha_x_temp)^2
      m_temp6              <- do.call(rbind, q_sig2_theta_is[[iter + 1]])
      sig2_theta_i[[iter + 1]] <- rowSums(m_temp5 + m_temp6) / (S + 2 * lambda)
      sig2_theta_i[[iter + 1]][in_sig2_theta] <- 0
      
      # Update sig2_beta
      sig2_beta[[iter + 1]] <- mean(q_mu_beta_0j[[iter + 1]]^2 +
                                      q_sig2_beta_0j[[iter + 1]])
      
      # Update sig2_b_j (with log penalty)
      m_temp8              <- (q_mu_b_js[[iter + 1]] - q_mu_beta_0j[[iter + 1]])^2 +
        q_sig2_beta_0j[[iter + 1]] + q_sig2_b_js[[iter + 1]]
      sig2_b_j[[iter + 1]] <- rowSums(m_temp8) / (S + 2 * lambda)
      sig2_b_j[[iter + 1]][in_sig2_b0] <- 0
      
      tau_M <- max(
        abs(alpha[[iter + 1]]        - alpha_old),
        abs(sig2_theta_i[[iter + 1]] - sig2_theta_i_old),
        abs(sig2_b_j[[iter + 1]]     - sig2_b_j_old),
        abs(sig2_beta[[iter + 1]]    - sig2_beta_old)
      )
    }
    
    # -- Convergence check ---------------------------------------------------
    tau <- max(
      abs(alpha[[iter + 1]]        - alpha[[iter]]),
      abs(sig2_theta_i[[iter + 1]] - sig2_theta_i[[iter]]),
      abs(sig2_b_j[[iter + 1]]     - sig2_b_j[[iter]]),
      abs(sig2_beta[[iter + 1]]    - sig2_beta[[iter]])
    )
    
    iter <- iter + 1
  }
  
  # ---------------------------------------------------------------------------
  # Post-convergence: sparsity thresholding
  # ---------------------------------------------------------------------------
  if (!debias) {
    sig2_b_j[[iter]][sig2_b_j[[iter]] < rho_N] <- 0
    out_sig2_b0    <- which(sig2_b_j[[iter]] != 0)
    in_sig2_b0     <- which(sig2_b_j[[iter]] == 0)
    
    sig2_theta_i[[iter]][sig2_theta_i[[iter]] < rho_N] <- 0
    out_sig2_theta <- which(sig2_theta_i[[iter]] != 0)
    in_sig2_theta  <- which(sig2_theta_i[[iter]] == 0)
  }
  
  # Apply variance floor (rho_N2) for likelihood computation
  sig2_b_j2     <- sig2_b_j[[iter]]
  sig2_theta_i2 <- sig2_theta_i[[iter]]
  
  sig2_b_j2[which((sig2_b_j2 < rho_N2) & (sig2_b_j2 != 0))]         <- rho_N2
  sig2_theta_i2[which((sig2_theta_i2 < rho_N2) & (sig2_theta_i2 != 0))] <- rho_N2
  
  # ---------------------------------------------------------------------------
  # Compute variational lower bound (likelihood) and model selection criteria
  # ---------------------------------------------------------------------------
  temp1_like    <- list()
  b_sum         <- q_mu_b_js[[iter]]^2 + q_sig2_b_js[[iter]]
  part7_temp    <- list()
  alpha_x_like  <- alpha[[iter]] %*% t(dummy_matrix_tilde)
  
  for (i in 1:I) {
    part1 <- -log(1 + exp(-sqrt(ksi_sij2[[iter]][[i]])))
    
    part2 <- (Ysij[[i]] - 1 / 2) *
      (t(-q_mu_b_js[[iter]]) + q_mu_theta_is[[iter]][[i]]) -
      (1 / 2) * sqrt(ksi_sij2[[iter]][[i]])
    
    part31 <- q_mu_theta_is[[iter]][[i]]^2 + q_sig2_theta_is[[iter]][[i]]
    part32 <- -2 * t(q_mu_b_js[[iter]]) * q_mu_theta_is[[iter]][[i]]
    part33 <- t(b_sum) - ksi_sij2[[iter]][[i]]
    part3  <- eta(ksi_sij2[[iter]][[i]]) * (part31 + part32 + part33)
    
    temp1_like[[i]] <- sum(part1 + part2 - part3)
    
    part7_temp[[i]] <- log(sig2_theta_i2[i]) +
      ((q_mu_theta_is[[iter]][[i]] - alpha_x_like[i, ])^2 +
         q_sig2_theta_is[[iter]][[i]]) / sig2_theta_i2[i]
  }
  
  part4 <- (S * (length(out_sig2_b0)) + J + S * (length(out_sig2_theta))) * log(2 * pi)
  
  part5 <- sum(log(sig2_beta[[iter]]) +
                 (q_mu_beta_0j[[iter]]^2 + q_sig2_beta_0j[[iter]]) / sig2_beta[[iter]])
  
  part6 <- sum(
    log(sig2_b_j2[out_sig2_b0]) +
      ((q_mu_b_js[[iter]][out_sig2_b0, ] - q_mu_beta_0j[[iter]][out_sig2_b0])^2 +
         q_sig2_b_js[[iter]][out_sig2_b0, ] +
         q_sig2_beta_0j[[iter]][out_sig2_b0]) / sig2_b_j2[out_sig2_b0]
  )
  
  part7 <- sum(do.call(rbind, part7_temp[out_sig2_theta]))
  
  temp2_like <- part4 + part5 + part6 + part7
  
  likelihood <- Reduce("+", temp1_like) - (1 / 2) * Reduce("+", temp2_like)
  
  N_total <- sum(Ns)
  df      <- length(out_sig2_b0) + length(out_sig2_theta)
  
  BIC <- -2 * (Reduce("+", temp1_like) - (1 / 2) * temp2_like) +
    df * log(N_total)
  
  GIC <- -2 * (Reduce("+", temp1_like) - (1 / 2) * temp2_like) +
    df * c * log(N_total) * log(log(N_total))
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  all_return <- list(
    sig2_beta    = sig2_beta[[iter]],
    sig2_theta_i = sig2_theta_i[[iter]],
    sig2_b_j     = sig2_b_j[[iter]],
    alpha        = alpha[[iter]],
    alpha_main   = alpha[[iter]][, -1],
    theta_i      = alpha[[iter]][, 1],
    likelihood   = likelihood,
    BIC          = BIC,
    GIC          = GIC,
    iter         = iter,
    tau          = tau
  )
  
  return(all_return)
}