Port_adjacency <- apply(as.matrix(read.csv(file = "adj_mat.csv")),2, as.numeric)
A_norm <- col_normalise_adjacency_matrix(Port_adjacency)
A_norm

burnin <- 100
N <- 1000
Ntotal <- N+burnin 
Delta <- 1
d <- 24

rate <- 5
c <- -0.8

# Choose the K matrix
identity_matrix <- diag(d)
K <- -identity_matrix - c * t(A_norm)

# Distributions
alpha_g <- 1.5
Q_dist <- function() stats::rgamma(1, shape = alpha_g, rate = 1) # Random 1-dimensional vector

U_dist <- function(d) {  rnorm(d) }


MCiterations <- 1000 

start_time <- proc.time()
for(i in 1:MCiterations){
  
  set.seed(1000+i)
  
  result <- simulate_graphsupOU(rate, K, Q_dist, U_dist, Delta, Ntotal, d)
  X <- result$X_t[(burnin+1):Ntotal,]
  
  file_name <- paste("POR_sim_data_", i, ".txt", sep = "")
  
  write.table(X, file = file_name, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  print(paste("Data written to:", file_name))
  if (i == 1) {
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    print(paste("Time for first iteration:", elapsed_time["elapsed"], "seconds"))
  }
  
}

