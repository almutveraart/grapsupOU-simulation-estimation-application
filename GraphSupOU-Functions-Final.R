library(expm) 
library(igraph) 
library(MASS) 

####Simulation
simulate_graphsupOU <- function(rate, K, Q_dist, U_dist, Delta, N, d) {
  
  T <- N * Delta
  
  # Simulate the Poisson arrival times
  arrival_times <- cumsum(rexp(ceiling(rate * T * 1.5), rate = rate)) 
  arrival_times <- arrival_times[arrival_times <= T] 
  
  # Simulate random memory parameter and jump size vector
  n_jumps <- length(arrival_times)
  Q_i <- replicate(n_jumps, Q_dist(), simplify = FALSE) 
  U_i <- lapply(1:n_jumps, function(i) U_dist(d)) 
  
  # Initialise results
  t_vals <- seq(Delta, N * Delta, by = Delta) 
  X_t <- matrix(0, nrow = length(t_vals), ncol = d)
  
  # Compute X(t) 
  for (j in seq_along(t_vals)) {
    t <- t_vals[j]
    relevant_indices <- which(arrival_times <= t) 
    tau <- arrival_times[relevant_indices]
    
    for (k in seq_along(relevant_indices)) {
      i <- relevant_indices[k]
      Q <- Q_i[[i]]
      U <- U_i[[i]]
      decay_time <- t - tau[k]
      X_t[j, ] <- X_t[j, ] + as.numeric(expm(Q * K* decay_time) %*% U)
    }
  }
  
  # Return time points and results
  return(list(time = t_vals, X_t = X_t))
}


#Normalise the adjacency matrix
col_normalise_adjacency_matrix <- function(A) {
  
  # Determine the number of nodes (columns of A)
  d <- ncol(A) 
  
  # Compute in-degrees n_j for each column
  in_degrees <- apply(A, 2, sum) 
  in_degrees <- pmax(1, in_degrees) 
  
  # Compute normalisation factors (1 / n_j)
  norm_factors <- 1 / in_degrees
  
  # Normalise columns
  A_normalised <- A %*% diag(norm_factors)
  
  return(A_normalised)
}



######Estimation
#Empirical covariance est:
autocovariance <- function(X, lag) {
  # X: d-dimensional time series matrix (N x d)
  # h: lag for autocovariance computation
  # on th egrid of the time series, i.e if the
  # time series is on a Delta grid,
  # we compute Cov(X_0,X_{h Delta})
  h<-lag
  
  N <- nrow(X)  # Number of observations (rows)
  
  # Compute the sample mean of X
  X_bar <- colMeans(X)
  
  # Initialize the autocovariance matrix
  autocovariance <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  
  # Compute the sum for the empirical autocovariance
  for (t in 1:(N - h)) {
    Xt <- X[t, , drop = FALSE]  
    Xt_h <- X[t + h, , drop = FALSE]  
    
    # Update the autocovariance sum
    autocovariance <- autocovariance + t(Xt-X_bar) %*% (Xt_h-X_bar) #note that Xt is a row and vector and needs to be transposed
  }
  
  # Compute the mean term
  #mean_term <- X_bar %*% t(X_bar)
  
  # Calculate the final autocovariance matrix
  autocovariance <- (1 / N) * autocovariance #- mean_term
  
  return(autocovariance)
}



# Function to compute maximum eigenvalues for lagged covariance matrices
find_max_eigen <- function(X, max_lag) {
  T <- nrow(X)  # Number of rows (time points)
  d <- ncol(X)  # Number of columns (time series dimension)
  
  # Compute Cov(X_t), which is the covariance of X_t with itself (lag 0)
  cov_X <- cov(X)
  inv_cov_X <- solve(cov_X) # this is cov_X^{-1}
  
  # Initialize vector to store max eigenvalues
  max_eigenvalues <- numeric(max_lag)
  max_real_eigenvalues <- numeric(max_lag)
  
  for (h in 1:max_lag) {
    
    # Compute autocovariance matrix for lag h
    autocov_h<-autocovariance(X, h)
    
    # Compute R(h; xi)
    R_h <- autocov_h %*% inv_cov_X
    
    
    # Compute eigenvalues and filter real eigenvalues
    eigenvalues <- eigen(R_h)$values
    real_eigenv <- Re(eigenvalues[Im(eigenvalues) == 0])
    
    # Take the eigenvalue with maximum real part
    max_eigenvalues[h] <- eigenvalues[which.max(Re(eigenvalues))]
    #
    if (length(real_eigenv) > 0) {
      max_real <- max(real_eigenv)
      max_real_eigenvalues[h]<- max_real
    } else {
      print("No real eigenvalue.")
    }
    
  }
  
  #return(max_eigenvalues) 
  return(list(max_eigenvalues=max_eigenvalues, max_real_eigenvalues=max_real_eigenvalues))
}


# Define objective function
objective_function <- function(params, h,  Delta, observed) {
  lambda <- params[1]
  alpha <- params[2]
  
  predicted <- (1+(1+lambda) * h *Delta)^(1-alpha)
  return(sum((observed - predicted)^2))  # Residual sum of squares
}




# Define the function
f <- function(h, lambda, alpha, Delta) {
  (1 +(1+lambda) * h * Delta)^(1 - alpha)
}


###OU case
# Define objective function
objective_function_OU <- function(params, h,  Delta, observed) {
  lambda <- params[1]
  
    predicted <- exp(-lambda* h *Delta)
  return(sum((observed - predicted)^2))  # Residual sum of squares
}

# Define the function
f_exp <- function(h, lambda, Delta) {
  exp(-lambda * h * Delta)
}

###OU case
# Define objective function
objective_function_2OU <- function(params, h,  Delta, observed) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  w <- params[3]
  
  predicted <- w*exp(-lambda1* h *Delta)+(1-w)*exp(-lambda2* h *Delta)
  return(sum((observed - predicted)^2))  # Residual sum of squares
}

# Define the function
f_2exp <- function(h, lambda1, lambda2, w, Delta) {
  w*exp(-lambda1* h *Delta)+(1-w)*exp(-lambda2* h *Delta)
}




#Compute K(c) as a function of Gamma
# Function to compute B
compute_K <- function(Gamma, h, alpha) {
  if (h == 0) stop("h cannot be zero.")
  
  # Compute the matrix power Gamma^(1/(1-alpha))
  exponent <- 1 / (1 - alpha)
  Gamma_power <- expm(logm(Gamma) * exponent)
  
  # Solve for K
  K <- (diag(nrow(Gamma)) - Gamma_power) / h
  
  return(K)
}





#Create K matrix
create_K <-function(c, A_norm){
  d <- nrow(A_norm)
  identity_matrix <- diag(d)
  K <- -identity_matrix - c * t(A_norm)
  return(K)
}

###############error measures
##Root mean square errors
rmse <- function(x, a){
  if (length(a) == 1) {
    a <- rep(a, length(x))  
  }
  sqrt(mean((x - a)^2))
}

median_absolute_error <- function(x, a) {
  if (length(a) == 1) {
    a <- rep(a, length(x))  
  }
  
  abs_error <- abs(a - x)
  
  return(median(abs_error))
}
