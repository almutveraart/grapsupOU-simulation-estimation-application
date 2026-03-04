# ============================================================
# GraphSupOU-Functions-Final.R
#
# Core functions for simulation, normalisation, estimation,
# and model comparison for the graph supOU process, as defined
# in the accompanying paper:
#
#   Mehta, S. and Veraart, A. E. D. (2025).
#   "Statistical inference for Lévy-driven graph supOU processes:
#    From short- to long-memory in high-dimensional time series."
#   arXiv:2502.08838 [stat.ME] (preprint).
#   https://arxiv.org/abs/2502.08838
#
# All equation references (e.g. "equation (4)") refer to the
# above preprint.
#
# Sections in this file:
#   1. Simulation
#   2. Graph normalisation
#   3. Estimation — empirical autocovariance and eigenvalues
#   4. Objective functions and fitted eigenvalue curves
#   5. Unconstrained reparametrisation (Gamma supOU model)
#   6. Auxiliary matrix functions
#   7. Error measures
#
# Dependencies: expm, igraph, MASS
# Minimum requirements: R >= 4.1.0
# ============================================================
library(expm)
library(igraph)
library(MASS)


# ============================================================
# 1. SIMULATION
# ============================================================
# simulate_graphsupOU()
#
# Simulates a discretely observed graph supOU process
#   X(t) = sum_{i: tau_i <= t} exp(Q_i * K * (t - tau_i)) * U_i
# at equally spaced times t = Delta, 2*Delta, ..., N*Delta.
#
# ACCELERATED VERSION: replaces repeated expm() calls with a one-time
# eigendecomposition of K. Since K is fixed across all jumps, we diagonalise
# it once:
#   K = V %*% diag(eig_vals) %*% V_inv
# Each matrix exponential is then evaluated as:
#   expm(Q_i * K * s) = V %*% diag(exp(Q_i * eig_vals * s)) %*% V_inv
# replacing an O(d^3) expm() call with O(d) scalar exponentials plus an
# O(d^2) matrix-vector product. Total cost: O(d^3 + varrho * N^2 * d^2),
# a factor of O(d) improvement over the original. For d=24, varrho=5, N=1000
# this amounts to a speedup of approximately 24 in floating-point operations.
# See equation (eq:fastsim) in the paper for the mathematical derivation.
#
# Arguments:
#   rate    -- Poisson arrival rate varrho of the background driving process
#   K       -- d x d drift matrix; typically K = create_K(c, A_norm)
#   Q_dist  -- zero-argument function returning a random scalar memory
#               parameter Q_i ~ pi (e.g. function() rgamma(1, alpha, 1))
#   U_dist  -- function(d) returning a random d-vector jump size U_i
#               (e.g. function(d) rnorm(d))
#   Delta   -- observation interval (time step)
#   N       -- number of observed time points; observations at Delta, ..., N*Delta
#   d       -- dimension of the process (number of nodes in the graph)
#
# Returns:
#   A list with components:
#     $time  -- length-N numeric vector of observation times (Delta, ..., N*Delta)
#     $X_t   -- N x d matrix; row j is the simulated process at time j*Delta
simulate_graphsupOU <- function(rate, K, Q_dist, U_dist, Delta, N, d) {
  
  T <- N * Delta
  
  # Simulate Poisson arrival times
  arrival_times <- cumsum(rexp(ceiling(rate * T * 1.5), rate = rate))
  arrival_times <- arrival_times[arrival_times <= T]
  
  # Simulate random memory parameters (scalars) and jump size vectors
  n_jumps <- length(arrival_times)
  Q_i     <- replicate(n_jumps, Q_dist(), simplify = TRUE)  # length-n_jumps numeric vector
  U_i     <- matrix(sapply(1:n_jumps, function(i) U_dist(d)), nrow = d)  # d x n_jumps matrix
  
  # Initialise output
  t_vals <- seq(Delta, N * Delta, by = Delta)
  X_t    <- matrix(0, nrow = length(t_vals), ncol = d)
  
  # Diagonalise K once: K = V %*% diag(eig_vals) %*% V_inv
  # eig_vals and columns of V may be complex if K is not symmetric
  eig      <- eigen(K)
  V        <- eig$vectors   # d x d eigenvector matrix
  eig_vals <- eig$values    # length-d eigenvalue vector (possibly complex)
  V_inv    <- solve(V)
  
  # Iterate over jumps (outer loop) and future time points (inner loop).
  # For each jump i, precompute V_inv %*% U_i[,i] (the jump vector in the
  # eigenbasis) once, then reuse it across all future time points.
  for (k in 1:n_jumps) {
    
    tau        <- arrival_times[k]
    future_idx <- which(t_vals >= tau)
    if (length(future_idx) == 0) next
    decay_times <- t_vals[future_idx] - tau
    
    # Rotate jump vector into eigenbasis once per jump: U_tilde = V^{-1} U_i
    V_inv_U <- V_inv %*% U_i[, k]
    
    for (jj in seq_along(future_idx)) {
      # e_ij = (exp(Q_i * k_1 * s), ..., exp(Q_i * k_d * s))  [elementwise]
      exp_diag <- exp(Q_i[k] * eig_vals * decay_times[jj])
      
      # Accumulate: X_{t_j} += Re(V %*% (e_ij * U_tilde_i))
      # Re() discards negligible imaginary parts from floating-point rounding;
      # the true result is real since expm() of a real matrix is real.
      X_t[future_idx[jj], ] <- X_t[future_idx[jj], ] +
        Re(as.numeric(V %*% (exp_diag * V_inv_U)))
    }
  }
  
  return(list(time = t_vals, X_t = X_t))
}


# ============================================================
# 2. GRAPH NORMALISATION
# ============================================================

# col_normalise_adjacency_matrix()
#
# Computes the column-normalised adjacency matrix
#   A_bar = A %*% D,
# where D = diag(1/n_1, ..., 1/n_d) and n_j = max(1, col_sum_j(A)).
# Isolated nodes (zero in-degree) are assigned n_j = 1 to avoid division
# by zero, leaving their column unchanged.
# See Section 2.3 of the paper.
#
# Arguments:
#   A -- d x d adjacency matrix (non-negative entries; need not be binary)
#
# Returns:
#   A d x d column-normalised matrix A_bar such that each column sums to
#   at most 1 (exactly 1 for non-isolated nodes).
col_normalise_adjacency_matrix <- function(A) {
  
  d            <- ncol(A)
  in_degrees   <- apply(A, 2, sum)
  in_degrees   <- pmax(1, in_degrees)   # guard against isolated nodes
  norm_factors <- 1 / in_degrees
  A_normalised <- A %*% diag(norm_factors)
  
  return(A_normalised)
}



# ============================================================
# 3. ESTIMATION — EMPIRICAL AUTOCOVARIANCE AND EIGENVALUES
# ============================================================

# autocovariance()
#
# Computes the empirical autocovariance matrix at lag h:
#   Gamma_hat(h) = (1 / (N - h)) * sum_{t=1}^{N-h} (X_t - X_bar)(X_{t+h} - X_bar)^T
#
# The divisor (N - h) matches the estimator in the paper (equation before (12))
# and differs from the R default cov() divisor of (N - 1).
#
# Arguments:
#   X   -- N x d matrix of observations (rows = time, columns = nodes)
#   lag -- integer lag h >= 1
#
# Returns:
#   A d x d symmetric matrix Gamma_hat(h).
autocovariance <- function(X, lag) {
  # X:   N x d matrix of observations
  # lag: integer lag h >= 1
  h <- lag
  N <- nrow(X)
  
  X_bar       <- colMeans(X)
  autocov_mat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  
  for (t in 1:(N - h)) {
    Xt   <- X[t,     , drop = FALSE]
    Xt_h <- X[t + h, , drop = FALSE]
    autocov_mat <- autocov_mat + t(Xt - X_bar) %*% (Xt_h - X_bar)
  }
  
  # Divide by (N - h) to match the paper's estimator
  autocov_mat <- (1 / (N - h)) * autocov_mat
  
  return(autocov_mat)
}


# find_max_eigen()
#
# Computes the largest real eigenvalue of the scaled empirical autocovariance
# matrix R_hat(h) = Gamma_hat(h) %*% Gamma_hat(0)^{-1} for lags h = 1, ..., max_lag.
#
# Two eigenvalue sequences are returned:
#   max_eigenvalues      -- eigenvalue with largest real part at each lag
#                           (may be complex due to numerical noise)
#   max_real_eigenvalues -- largest purely real eigenvalue at each lag;
#                           this is the quantity used for estimation (equation (12)
#                           and surrounding discussion in Section 3.1 of the paper)
#
# At lags where the empirical autocovariance estimate is noisy and no purely
# real eigenvalue exists, max_real_eigenvalues[h] is set to NA. These entries
# are silently dropped by objective_function() so they do not influence the
# parameter estimates.
#
# Arguments:
#   X       -- N x d matrix of observations
#   max_lag -- maximum lag H; eigenvalues computed for h = 1, ..., max_lag
#
# Returns:
#   A list with components:
#     $max_eigenvalues      -- length-max_lag numeric vector (possibly complex)
#     $max_real_eigenvalues -- length-max_lag numeric vector (NA where unavailable)
find_max_eigen <- function(X, max_lag) {
  
  N <- nrow(X)
  d <- ncol(X)
  
  # Lag-0 covariance and its inverse (used for all lags)
  cov_X     <- cov(X)
  inv_cov_X <- solve(cov_X)
  
  max_eigenvalues      <- numeric(max_lag)
  max_real_eigenvalues <- rep(NA_real_, max_lag)  # NA until filled; avoids spurious 0s
  
  for (h in 1:max_lag) {
    
    # Scaled empirical autocovariance: R_hat(h) = cov_hat(h) %*% cov_hat(0)^{-1}
    autocov_h <- autocovariance(X, h)
    R_h       <- autocov_h %*% inv_cov_X
    
    # Eigenvalues of R_hat(h)
    eigenvalues <- eigen(R_h)$values
    real_eigenv <- Re(eigenvalues[Im(eigenvalues) == 0])
    
    # Eigenvalue with largest real part (may be complex)
    max_eigenvalues[h] <- eigenvalues[which.max(Re(eigenvalues))]
    
    # Largest purely real eigenvalue (used for estimation).
    # If no purely real eigenvalue exists (can occur at large lags due to noise
    # in the autocovariance estimate), store NA. The objective_function drops
    # NA entries so they do not influence the parameter estimates.
    if (length(real_eigenv) > 0) {
      max_real_eigenvalues[h] <- max(real_eigenv)
    } else {
      max_real_eigenvalues[h] <- NA_real_
    }
  }
  
  return(list(max_eigenvalues      = max_eigenvalues,
              max_real_eigenvalues = max_real_eigenvalues))
}



# ============================================================
# 4. OBJECTIVE FUNCTIONS AND FITTED EIGENVALUE CURVES
# ============================================================

# objective_function()
#
# Mean squared error loss for the Gamma supOU model (equation (11) in the paper):
#   rho(h; alpha, c) = (1 + (1 + c) * h * Delta)^(1 - alpha)
#
# Lags where observed is NA (no real eigenvalue found) are silently dropped.
#
# Arguments:
#   params   -- numeric vector of length 2: params[1] = c (network parameter,
#               equals c * a* with a* = 1 for connected graphs; named lambda
#               internally for historical reasons), params[2] = alpha (Gamma
#               shape parameter, must satisfy alpha > 1)
#   h        -- integer vector of lags h = 1, ..., H
#   Delta    -- observation interval
#   observed -- numeric vector of observed max real eigenvalues (from find_max_eigen)
#
# Returns:
#   Scalar MSE between observed and model-predicted eigenvalue curves.
objective_function <- function(params, h, Delta, observed) {
  lambda <- params[1]   # = c (named lambda for historical reasons)
  alpha  <- params[2]
  # Drop lags where no real eigenvalue was found (stored as NA)
  valid    <- !is.na(observed)
  h        <- h[valid]
  observed <- observed[valid]
  predicted <- (1 + (1 + lambda) * h * Delta)^(1 - alpha)
  return(mean((observed - predicted)^2))
}


# f()
#
# Fitted eigenvalue curve for the Gamma supOU model; vectorised over h.
# Companion to objective_function() for plotting fitted curves.
#
# Arguments:
#   h      -- numeric vector of lags
#   lambda -- estimated network parameter c (see objective_function)
#   alpha  -- estimated Gamma shape parameter
#   Delta  -- observation interval
#
# Returns:
#   Numeric vector of predicted eigenvalues (1 + (1 + lambda) * h * Delta)^(1 - alpha).
f <- function(h, lambda, alpha, Delta) {
  (1 + (1 + lambda) * h * Delta)^(1 - alpha)
}


# objective_function_OU()
#
# Mean squared error loss for the graph OU model (exponentially decaying
# eigenvalue function):
#   rho(h) = exp(-lambda * h * Delta)
# where lambda absorbs the factor (1 + c * a*) into a single decay rate.
#
# Arguments:
#   params   -- numeric vector of length 1: params[1] = lambda (decay rate > 0)
#   h        -- integer vector of lags
#   Delta    -- observation interval
#   observed -- numeric vector of observed max real eigenvalues
#
# Returns:
#   Scalar MSE between observed and model-predicted eigenvalue curves.
objective_function_OU <- function(params, h, Delta, observed) {
  lambda    <- params[1]
  predicted <- exp(-lambda * h * Delta)
  return(mean((observed - predicted)^2))
}


# f_exp()
#
# Fitted eigenvalue curve for the graph OU model; vectorised over h.
# Companion to objective_function_OU() for plotting fitted curves.
#
# Arguments:
#   h      -- numeric vector of lags
#   lambda -- estimated decay rate
#   Delta  -- observation interval
#
# Returns:
#   Numeric vector exp(-lambda * h * Delta).
f_exp <- function(h, lambda, Delta) {
  exp(-lambda * h * Delta)
}


# objective_function_2OU()
#
# Mean squared error loss for the sum-of-two-exponentials supOU model
# (equation (10) in the paper, two terms):
#   rho(h) = w * exp(-lambda1 * h * Delta) + (1 - w) * exp(-lambda2 * h * Delta)
#
# Arguments:
#   params   -- numeric vector of length 3:
#                 params[1] = lambda1 (first decay rate)
#                 params[2] = lambda2 (second decay rate)
#                 params[3] = w       (mixture weight in [0, 1])
#   h        -- integer vector of lags
#   Delta    -- observation interval
#   observed -- numeric vector of observed max real eigenvalues
#
# Returns:
#   Scalar MSE between observed and model-predicted eigenvalue curves.
objective_function_2OU <- function(params, h, Delta, observed) {
  lambda1   <- params[1]
  lambda2   <- params[2]
  w         <- params[3]
  predicted <- w * exp(-lambda1 * h * Delta) + (1 - w) * exp(-lambda2 * h * Delta)
  return(mean((observed - predicted)^2))
}


# f_2exp()
#
# Fitted eigenvalue curve for the sum-of-two-exponentials supOU model;
# vectorised over h. Companion to objective_function_2OU() for plotting.
#
# Arguments:
#   h       -- numeric vector of lags
#   lambda1 -- first decay rate
#   lambda2 -- second decay rate
#   w       -- mixture weight in [0, 1]
#   Delta   -- observation interval
#
# Returns:
#   Numeric vector w * exp(-lambda1*h*Delta) + (1-w) * exp(-lambda2*h*Delta).
f_2exp <- function(h, lambda1, lambda2, w, Delta) {
  w * exp(-lambda1 * h * Delta) + (1 - w) * exp(-lambda2 * h * Delta)
}


# ============================================================
# 5. UNCONSTRAINED REPARAMETRISATION (GAMMA supOU MODEL)
# ============================================================
#
# The natural parameter space has two constraints:
#   alpha > 1      and      c in (-1, 1)
#
# We remove both via smooth bijections:
#   phi = log(alpha - 1)    in R    [inverse: alpha = 1 + exp(phi)]
#   psi = arctanh(c)        in R    [inverse: c = tanh(psi)       ]
#
# Working in (phi, psi) space:
#   - eliminates all bound constraints
#   - removes the boundary pile-up near c -> -1 visible in contour plots
#   - allows gradient-based unconstrained optimisers (BFGS, Brent)
#
# After optim(), back-transform via backtransform_params() to recover
# (c, alpha) on their natural scales.


# transform_params()
#
# Forward transformation: (c, alpha) --> (psi, phi).
#
# Arguments:
#   c_val     -- network parameter c in (-1, 1)
#   alpha_val -- Gamma shape parameter alpha > 1
#
# Returns:
#   Numeric vector c(psi, phi) on the unconstrained real scale.
transform_params <- function(c_val, alpha_val) {
  psi <- atanh(c_val)                # c in (-1,1) --> psi in R
  phi <- log(alpha_val - 1)          # alpha > 1   --> phi in R
  c(psi, phi)
}


# backtransform_params()
#
# Inverse transformation: (psi, phi) --> (c, alpha).
#
# Arguments:
#   psi -- unconstrained transform of c (real scalar)
#   phi -- unconstrained transform of alpha (real scalar)
#
# Returns:
#   A list with components $c (in (-1,1)) and $alpha (> 1).
backtransform_params <- function(psi, phi) {
  c_val     <- tanh(psi)             # psi in R --> c in (-1,1)
  alpha_val <- 1 + exp(phi)          # phi in R --> alpha > 1
  list(c = c_val, alpha = alpha_val)
}


# objective_function_unconstrained()
#
# Unconstrained objective function for joint optimisation over (psi, phi).
# Wraps objective_function() by back-transforming (psi, phi) to (c, alpha)
# before evaluating the MSE.
#
# Arguments:
#   params   -- numeric vector of length 2: params[1] = psi, params[2] = phi
#   h        -- integer vector of lags
#   Delta    -- observation interval
#   observed -- numeric vector of observed max real eigenvalues
#
# Returns:
#   Scalar MSE (identical to objective_function() evaluated at back-transformed params).
objective_function_unconstrained <- function(params, h, Delta, observed) {
  psi <- params[1]
  phi <- params[2]
  bt  <- backtransform_params(psi, phi)
  # Reuse the constrained objective on the back-transformed values
  objective_function(c(bt$c, bt$alpha), h = h, Delta = Delta, observed = observed)
}


# objective_function_c_unconstrained()
#
# Unconstrained objective function for Step 2 of the two-step procedure:
# only psi (transformed c) is free; phi (transformed alpha) is fixed at
# phi_fixed = log(alpha_hat - 1) from Step 1.
#
# Arguments:
#   psi       -- unconstrained transform of c (optimised by optim)
#   phi_fixed -- fixed value of phi from Step 1
#   h         -- integer vector of lags
#   Delta     -- observation interval
#   observed  -- numeric vector of observed max real eigenvalues
#
# Returns:
#   Scalar MSE at the given psi with alpha held fixed.
objective_function_c_unconstrained <- function(psi, phi_fixed, h, Delta, observed) {
  bt <- backtransform_params(psi, phi_fixed)
  objective_function(c(bt$c, bt$alpha), h = h, Delta = Delta, observed = observed)
}


# fit_graphsupOU_unconstrained()
#
# Convenience wrapper for joint unconstrained optimisation over (psi, phi).
# Runs optim() with the specified method (default: BFGS) and returns
# estimates on the original (c, alpha) scale together with convergence info.
#
# Arguments:
#   observed -- numeric vector of observed max real eigenvalues
#   h        -- integer vector of lags
#   Delta    -- observation interval
#   c0       -- starting value for c (default 0.5); transformed to psi0 = atanh(c0)
#   alpha0   -- starting value for alpha (default 2); transformed to phi0 = log(alpha0-1)
#   method   -- optim() method (default "BFGS")
#
# Returns:
#   A list with components:
#     $c         -- estimated network parameter on (-1, 1)
#     $alpha     -- estimated Gamma shape parameter (> 1)
#     $psi_hat   -- estimated psi (unconstrained transform of c)
#     $phi_hat   -- estimated phi (unconstrained transform of alpha)
#     $value     -- minimised MSE
#     $converged -- logical; TRUE if optim() convergence code is 0
fit_graphsupOU_unconstrained <- function(observed, h, Delta,
                                         c0 = 0.5, alpha0 = 2,
                                         method = "BFGS") {
  psi0 <- atanh(c0)
  phi0 <- log(alpha0 - 1)
  
  fit <- optim(
    par    = c(psi0, phi0),
    fn     = objective_function_unconstrained,
    h      = h,
    Delta  = Delta,
    observed = observed,
    method = method,
    control  = list(maxit = 1000)
  )
  
  bt <- backtransform_params(fit$par[1], fit$par[2])
  list(
    c         = bt$c,
    alpha     = bt$alpha,
    psi_hat   = fit$par[1],
    phi_hat   = fit$par[2],
    value     = fit$value,
    converged = (fit$convergence == 0)
  )
}


# fit_c_unconstrained()
#
# Convenience wrapper for Step 2 of the two-step procedure: unconstrained
# optimisation over c only, with alpha fixed at alpha_fixed.
# Uses Brent's method (one-dimensional) over psi in [-6, 6], corresponding
# to c in (tanh(-6), tanh(6)) ≈ (-0.9999, 0.9999).
#
# Arguments:
#   observed    -- numeric vector of observed max real eigenvalues
#   h           -- integer vector of lags
#   Delta       -- observation interval
#   alpha_fixed -- fixed alpha estimate from Step 1 (must be > 1)
#   c0          -- starting value for c (default 0.5)
#
# Returns:
#   A list with components:
#     $c         -- estimated network parameter on (-1, 1)
#     $psi_hat   -- estimated psi (unconstrained transform of c)
#     $value     -- minimised MSE
#     $converged -- logical; TRUE if optim() convergence code is 0
fit_c_unconstrained <- function(observed, h, Delta, alpha_fixed, c0 = 0.5) {
  phi_fixed <- log(alpha_fixed - 1)
  psi0      <- atanh(c0)
  
  fit <- optim(
    par      = psi0,
    fn       = objective_function_c_unconstrained,
    phi_fixed = phi_fixed,
    h        = h,
    Delta    = Delta,
    observed = observed,
    method   = "Brent",
    lower    = -6,    # corresponds to c = tanh(-6) ≈ -0.9999
    upper    =  6     # corresponds to c = tanh( 6) ≈  0.9999
  )
  
  list(
    c         = tanh(fit$par),
    psi_hat   = fit$par,
    value     = fit$value,
    converged = (fit$convergence == 0)
  )
}

# ============================================================
# 6. AUXILIARY MATRIX FUNCTIONS
# ============================================================

# create_K()
#
# Constructs the drift matrix K(c) = -I - c * A_bar^T.
# See equation (4) in Section 2.4 of the paper.
#
# Arguments:
#   c      -- network parameter c in (-1, 1)
#   A_norm -- d x d column-normalised adjacency matrix (from col_normalise_adjacency_matrix)
#
# Returns:
#   d x d drift matrix K(c).
create_K <- function(c, A_norm) {
  d <- nrow(A_norm)
  K <- -diag(d) - c * t(A_norm)
  return(K)
}

# ============================================================
# 7. ERROR MEASURES
# ============================================================

# rmse()
#
# Root mean squared error between a vector of estimates x and true value(s) a.
#
# Arguments:
#   x -- numeric vector of estimates
#   a -- numeric scalar (recycled) or vector of true values
#
# Returns:
#   Scalar RMSE: sqrt(mean((x - a)^2)).
rmse <- function(x, a) {
  if (length(a) == 1) a <- rep(a, length(x))
  sqrt(mean((x - a)^2))
}


# median_absolute_error()
#
# Median absolute error between a vector of estimates x and true value(s) a.
#
# Arguments:
#   x -- numeric vector of estimates
#   a -- numeric scalar (recycled) or vector of true values
#
# Returns:
#   Scalar MAE: median(|x - a|).
median_absolute_error <- function(x, a) {
  if (length(a) == 1) a <- rep(a, length(x))
  median(abs(a - x))
}

# Median absolute error
median_absolute_error <- function(x, a) {
  if (length(a) == 1) a <- rep(a, length(x))
  median(abs(a - x))
}