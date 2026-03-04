# ============================================================
# Portugal-Sim.R
#
# Monte Carlo data generation for the simulation study in
# Section 4 of:
#
#   Mehta, S. and Veraart, A. E. D. (2026).
#   "Statistical inference for Lévy-driven graph supOU processes:
#    From short- to long-memory in high-dimensional time series."
#   arXiv:2502.08838 [stat.ME] (preprint).
#   https://arxiv.org/abs/2502.08838
#
# Simulates MCiterations = 1000 independent realisations of the
# graph supOU process (accelerated exact simulation, Section 4.1)
# on the 24-node Portuguese subnetwork, using the long-memory
# Gamma(alpha, 1) specification with parameters motivated by the
# empirical study (Section 5):
#   alpha = 1.5,  c = -0.8,  rate = 5,  Delta = 1,  N = 1000
#
# Each realisation is saved as a tab-separated text file in the
# simdata/ subdirectory. The analysis script Portugal-Sim-Analysis.R
# reads these files to produce the figures in Section 4.
#
# Required input files (must be in the working directory):
#   adj_mat.csv                  -- 24x24 adjacency matrix for
#                                   the Portuguese subnetwork
#   GraphSupOU-Functions-Final.R -- defines simulate_graphsupOU,
#                                   col_normalise_adjacency_matrix
#
# Output:
#   simdata/POR_sim_data_1.txt, ..., simdata/POR_sim_data_1000.txt
#   Each file: N x d tab-separated matrix (no row/col names),
#   rows = time points, columns = nodes.
#
# Parallelisation: uses all available cores minus one (via
#   parallel + doParallel + foreach).
#
# Reproducibility: each iteration uses set.seed(1000 + i),
#   ensuring results are fully reproducible across runs.
#
# Dependencies: parallel, doParallel, foreach, expm
# Minimum requirements: R >= 4.1.0
# ============================================================
source("GraphSupOU-Functions-Final.R")

library(parallel)
library(doParallel)
library(foreach)

# ── Graph structure ───────────────────────────────────────────────────────────
# Load the 24-node Portuguese subnetwork adjacency matrix and
# column-normalise it to obtain A_bar (Section 2.3, equation before (2))
Port_adjacency <- apply(as.matrix(read.csv(file = "adj_mat.csv")), 2, as.numeric)
A_norm <- col_normalise_adjacency_matrix(Port_adjacency)

# ── Simulation parameters ─────────────────────────────────────────────────────
# Parameters match those used in the empirical study (Section 5):
#   alpha = 1.5  --> long-memory regime (alpha in (1,2))
#   c     = -0.8 --> strong negative network effect
#   rate  = 5    --> Poisson arrival rate varrho of the compound Poisson basis
#   Delta = 1    --> observation interval (hourly, as in the empirical data)
#   N     = 1000 --> number of observations per simulation run
#   d     = 24   --> number of nodes
# A burn-in of 100 time steps is discarded to approximate stationarity.
burnin   <- 100
N        <- 1000
Ntotal   <- N + burnin
Delta    <- 1
d        <- 24
rate     <- 5
c        <- -0.8

# Drift matrix K(c) = -I - c * A_bar^T  (equation (5))
K <- -diag(d) - c * t(A_norm)

# ── Distributions for the driving Lévy basis ──────────────────────────────────
# Gamma(alpha, 1) memory distribution pi (equation (8)):
#   theta_2 ~ Gamma(alpha_g, 1),  alpha_g = 1.5
# Jump size distribution: U_i ~ N(0, I_d), giving mu_L = 0, sigma^2_L = 5*I_d
alpha_g <- 1.5
Q_dist  <- function() stats::rgamma(1, shape = alpha_g, rate = 1)
U_dist  <- function(d) rnorm(d)

MCiterations <- 1000

# ── Output directory ──────────────────────────────────────────────────────────
# Files are written to simdata/ as tab-separated N x d matrices.
# Portugal-Sim-Analysis.R expects files named POR_sim_data_i.txt.
dir.create("simdata", showWarnings = FALSE)

# ── Parallelisation setup ─────────────────────────────────────────────────────
# Detect available cores and leave one free for the OS
n_cores <- detectCores() - 1
cat(paste("Using", n_cores, "cores\n"))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export everything each worker needs
clusterExport(cl, varlist = c("simulate_graphsupOU", "K", "Q_dist", "U_dist",
                              "Delta", "Ntotal", "d", "rate", "burnin", "N"))

# Each worker also needs the expm library for simulate_graphsupOU
clusterEvalQ(cl, library(expm))

# ── Monte Carlo simulation loop ───────────────────────────────────────────────
# Each iteration:
#   1. Sets a deterministic seed (1000 + i) for reproducibility
#   2. Simulates Ntotal = N + burnin time steps
#   3. Discards the first burnin rows to approximate the stationary distribution
#   4. Writes the remaining N x d matrix to simdata/POR_sim_data_i.txt
start_time <- proc.time()

foreach(i = 1:MCiterations, .packages = "expm") %dopar% {
  
  # Set per-iteration seed for reproducibility
  set.seed(1000 + i)
  
  result <- simulate_graphsupOU(rate, K, Q_dist, U_dist, Delta, Ntotal, d)
  X      <- result$X_t[(burnin + 1):Ntotal, ]
  
  file_name <- paste("simdata/POR_sim_data_", i, ".txt", sep = "")
  write.table(X, file = file_name, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Return iteration number so foreach has something to collect
  i
}

end_time   <- proc.time()
elapsed    <- (end_time - start_time)["elapsed"]
cat(paste("Total time:", round(elapsed, 1), "seconds\n"))
cat(paste("Average time per iteration:", round(elapsed / MCiterations, 2), "seconds\n"))

# Shut down the cluster cleanly
stopCluster(cl)