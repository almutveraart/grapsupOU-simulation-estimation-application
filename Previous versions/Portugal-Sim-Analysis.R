library(ggplot2)
library(maps)
library(igraph)
library(tidyr)
library(reshape2)
library(dplyr)
library(patchwork)

#Read in graph structure
Port_adjacency <- apply(as.matrix(read.csv(file = "adj_mat.csv")),2, as.numeric)
A_norm <- col_normalise_adjacency_matrix(Port_adjacency)
A_norm


MonteCarloNumber <- 1000


#########################
#######################
# Estimate for different lags
# Create the list 
my_parameter_list <- vector("list", length = MonteCarloNumber)


for (i in 1:MonteCarloNumber) {
  
  file_name <- paste0("POR_sim_data_", i, ".txt")
  
  my_data <- read.table(file_name, header = FALSE, sep = "\t")
  
  X<-as.matrix(my_data)
  
  lag_vector <- seq(5, 100, by=1)
  l <- length(lag_vector)
  alpha_vector <- numeric(l)
  c_vector <- numeric(l)
  
  ev_all <- find_max_eigen(X, 100)
  max_real_eigenvalues_all <- ev_all$max_real_eigenvalues
  
  for (j in 1:l) {
    
    max_lag <- lag_vector[j]
    max_real_eigenvalues <- max_real_eigenvalues_all[1:max_lag]
    
    # Perform nonlinear least squares estimation
    Delta <- 1
    h_values <- (1:max_lag)*Delta
    initial_guess <- c(0.5, 2)  
    fit <- optim(
      par = initial_guess,
      fn = objective_function,
      h = h_values,
      Delta=Delta,
      observed = max_real_eigenvalues,
      method = "L-BFGS-B",
      lower = c(-0.9999, 1),  
      upper = c(0.9999, 10)  
    )
    
    c_vector[j] <- fit$par[1]
    alpha_vector[j] <- fit$par[2]
  }
  
  my_parameter_list[[i]] <- list(ev = max_real_eigenvalues_all[5:100], alpha_vector = alpha_vector, c_vector = c_vector)
}

my_parameter_list

##Convert to matrix
my_col <- 96
alpha_matrix <- matrix(0, nrow=MonteCarloNumber, ncol=my_col)
c_matrix <- matrix(0, nrow=MonteCarloNumber, ncol=my_col)

for(i in 1:MonteCarloNumber){
  alpha_matrix[i,] <- my_parameter_list[[i]]$alpha
  c_matrix[i,] <- my_parameter_list[[i]]$c
}

#Plot estimates in boxplots

alpha_df <- as.data.frame(alpha_matrix)
c_df <- as.data.frame(c_matrix)


alpha_long <- melt(alpha_df, variable.name = "Column", value.name = "Alpha")
c_long <- melt(c_df, variable.name = "Column", value.name = "C")


alpha_long$Column <- as.factor(alpha_long$Column)
c_long$Column <- as.factor(c_long$Column)


new_labels <- seq(5, 100)  
original_labels <- paste0("V", seq(1, 96))  
label_map <- setNames(new_labels, original_labels)  


custom_labels <- function(labels) {
  mapped_labels <- label_map[labels]  
  return(ifelse(mapped_labels %% 10 == 0, mapped_labels, ""))  
}


alpha_plot <- ggplot(alpha_long, aes(x = Column, y = Alpha)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA)+
  scale_x_discrete(labels = custom_labels) +
  labs(title = "", x = "Lag", y = expression(alpha)) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )+
  geom_hline(yintercept = 1.5,  color = "red")

c_plot <- ggplot(c_long, aes(x = Column, y = C)) +
  geom_boxplot(fill = "lightgreen", outlier.shape = NA)+
  scale_x_discrete(labels = custom_labels) +
  labs(title = "", x = "Lag", y = "c") +
   theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )+
  geom_hline(yintercept = -0.8,  color = "red")

combined_plot <- alpha_plot / c_plot
print(combined_plot)

ggsave("sim_alpha.eps",  plot=alpha_plot, width = 20, height = 20, units = "cm")
ggsave("sim_c.eps",  plot=c_plot, width = 20, height = 20, units = "cm")

################
#####
#Compute error measures

mean_alpha <- apply(alpha_matrix, 2, mean)
mean_alpha


mean_c <- apply(c_matrix, 2, mean)
mean_c



lags <- seq_along(mean_alpha)        # Lag indices (1 to 96 corresponding to 5 to 100)

data <- data.frame(
  lag = lags,
  med_alpha = apply(alpha_matrix, 2,median_absolute_error, a=1.5),
  rmse_alpha = apply(alpha_matrix, 2,rmse, a=1.5),
  med_c = apply(c_matrix, 2,median_absolute_error, a=-0.8),
  rmse_c = apply(c_matrix, 2,rmse, a=-0.8)
)

data_long <- data %>%
  pivot_longer(cols = c(med_alpha, rmse_alpha, med_c, rmse_c), 
               names_to = "variable", 
               values_to = "error")

p_error <-ggplot(data_long, aes(x = lag, y = error, color = variable)) +
   geom_line(data = subset(data_long, variable %in% c("med_alpha")), size = 1, linetype = "solid") +
  geom_line(data = subset(data_long, variable %in% c("med_c")), size = 1, linetype = "longdash") +
    geom_line(data = subset(data_long, variable %in% c("rmse_alpha")), size = 1, linetype = "dotted") +
  geom_line(data = subset(data_long, variable %in% c("rmse_c")), size = 1, linetype = "dotdash") +
  
  labs(
    title = "",
    x = "Lag",
    y = "Error measure",
    color = "Variable"
  ) +
  scale_color_manual(
    values = c("med_alpha" = "darkblue", "rmse_alpha"="blue", "med_c" = "darkgreen", "rmse_c"="green"),  # Set colors for lines
    labels = c(expression(alpha - median),  expression(c - median), expression(alpha - RMSE), expression(c - RMSE))  # Set legend labels with LaTeX-style alpha
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
p_error
ggsave("sim-error.eps",  plot=p_error, width = 20, height = 20, units = "cm")


#Provide mean and median for certain fixed lag
max_lag <- 35

alpha_vector <- alpha_matrix[, (max_lag-5+1)]
c_vector <- c_matrix[, (max_lag-5+1)]



data <- data.frame(
  value = c(alpha_vector, c_vector),
  group = rep(c("alpha", "c"), each = length(alpha_vector))
)

p_vio <- ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +               
  geom_hline(yintercept = 1.5,              
             color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -0.8,             
             color = "red", linetype = "dashed") +
  labs(title = "",
       x = "Parameter", y = "Value", fill="") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )+
  scale_x_discrete(
    labels = c(expression(alpha), "c")  
  ) +
  scale_fill_manual(values = c("alpha" = "skyblue", "c" = "pink"),
                    labels = c(expression(alpha), "c"))  

p_vio
ggsave("p_vio.eps", plot = p_vio, width = 20, height = 20, units = "cm")

#####Recover the mean and variance for each Monte Carlo run
my_mean_var_list <- vector("list", length = MonteCarloNumber)

average_mean <- numeric(24)
average_var <- matrix(0,24,24)

diag_var_all <- matrix(0,MonteCarloNumber,24)

for (i in 1:MonteCarloNumber) {
  
  print(i)
  
  file_name <- paste0("POR_sim_data_", i, ".txt")
  
  my_data <- read.table(file_name, header = FALSE, sep = "\t")
  
  X<-as.matrix(my_data)
  
  c <- c_vector[i]
  alpha <- alpha_vector[i]
  
  K <- create_K(c, A_norm)
  
  mean <- (1-alpha)*K%*%as.numeric(colMeans(X))
  var <- (1-alpha)*(K%*%cov(X)+cov(X)%*%t(K))
  
  diag_var_all[i, ] <- diag(var)
  
  my_mean_var_list[[i]] <- list(mean = mean, var = var)
  
  average_mean <- average_mean + mean/MonteCarloNumber
  average_var <- average_var + var/MonteCarloNumber
  
}

#Compute  the median values of the means and variance matrix elements


median_var_matrix <- matrix(0,24,24)
median_mu_vec <- numeric(24)

for(i in 1:24){
  for(j in 1:24){
    
    values <- sapply(1:1000, function(k) my_mean_var_list[[k]]$var[i, j])
    
    
    median_var_matrix[i, j] <- median(values)
  }
  mean_values <- sapply(1:1000, function(k) my_mean_var_list[[k]]$mean[i])
  median_mu_vec[i] <- median(mean_values)
}

var_data <- melt(median_var_matrix)
colnames(var_data) <- c("Row", "Column", "Value")

var_data$Row <- factor(var_data$Row)
var_data$Column <- factor(var_data$Column)

p<-ggplot(var_data, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.1, 5.1), name = "Cov") +
  labs(
    x = "Component",
    y = "Component"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )+
  scale_x_discrete(breaks = seq(4, 24, by = 4), labels = as.character(seq(4, 24, by = 4))) +
  scale_y_discrete(breaks = paste0("X", seq(4, 24, by = 4)), labels = seq(4, 24, by = 4))+
  geom_text(data = var_data %>% filter(Value > 4), aes(label = round(Value, 2)), size = 5, color = "black")



p
ggsave("estvar_median_sim.eps",  plot=p, width = 20, height = 20, units = "cm")
###################
#median of means plot
mean_data <- data.frame(Dimension = 1:24, Mean = median_mu_vec)

p <- ggplot(mean_data, aes(x = Dimension, y = Mean, fill = Mean)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.02, 0.02), name = "Mean") +
  labs(
    x = "Component",
    y = "Mean"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )

p

ggsave("estmean_median_sim.eps", plot = p, width = 20, height = 20, units = "cm")

