library(forecast)
library(ggplot2)
library(reshape2)
library(tidyr)

# Read in the data
data <- read.csv("wind_signal_COSMO.csv")

Time <- data[,1]
POR_data <- data[,2:25]

N <- nrow(POR_data) 

# View the first few rows
head(POR_data)

# We look at some exploratory plots
plot(POR_data$X1, type="l")
acf(ts(POR_data$X1))


# Prepare the data for plotting
Time_ad <- as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S")

POR_data_combined <- cbind(Time = Time_ad, POR_data)

POR_data_long <- melt(POR_data_combined, id.vars = "Time")

POR_data_long$variable <- as.numeric(gsub("X", "", POR_data_long$variable))

breaks <- seq(1, max(POR_data_long$variable), by = 4)

ggplot(POR_data_long, aes(x = variable, y = Time, fill = value)) +
  geom_tile() +
  scale_x_continuous(breaks = breaks) +  
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),    
    values = scales::rescale(c(-1, 0, 1)), 
    limits = c(-1, 1),                     
    name = "Value"
  ) +
  labs(
    x = "Component",
    y = "Time"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
ggsave("heatmap-original.eps", width = 20, height = 20, units = "cm")


##########################################
# Plot acfs of original data
lag_length <- 48

acf_list <- lapply(1:24, function(i) stats::acf(POR_data[, i], lag.max = lag_length, plot = FALSE))

acf_data <- data.frame(
  Lag = rep(0:lag_length, times = 24),  
  ACF_Value = unlist(lapply(acf_list, function(acf_result) acf_result$acf)),
  Column = rep(1:24, each = (lag_length+1))  
)

conf_limit <- 1.96 / sqrt(N)

ggplot(acf_data, aes(x = factor(Lag), y = ACF_Value)) +
  geom_boxplot() +
  labs(title = "", x = "Lag", y = "ACF") +
  scale_x_discrete(breaks = seq(0, lag_length, by = 24))+ 
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-conf_limit, conf_limit), 
             linetype = "dotted", color = "blue")+
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )


ggsave("empacfs-original.eps", width = 20, height = 20, units = "cm")

###################
# Plot the Lisbon data
Lisbon_data <- POR_data[, 22] 

time_series_data <- data.frame(Time_ad, Lisbon_data)

ggplot(time_series_data, aes(x = Time_ad, y = Lisbon_data)) +
  geom_line(color = "black", linewidth = 1) +
  labs( x = "Time", y = "Lisbon time series data") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
ggsave("Lisbondata-original.eps", width = 20, height = 20, units = "cm")



#####################################################
# Create a matrix to store the stationary components
stationary_matrix <- matrix(NA, nrow = nrow(POR_data), ncol = ncol(POR_data))

# Loop through each column of the dataset
for (i in seq_along(POR_data)) {
  
  ts_data <- POR_data[[i]]
  
  daily_period <- 24  
  yearly_period <- 24 * 365  
  
  # Decompose for daily seasonality
  stl_daily <- stl(ts(ts_data, frequency = daily_period), s.window = "periodic")
  deseasonalised_daily <- stl_daily$time.series[, "remainder"] + stl_daily$time.series[, "trend"]
  
  # Decompose for yearly seasonality
  stl_yearly <- stl(ts(deseasonalised_daily, frequency = yearly_period), s.window = "periodic")
  
  deseasonalised <- stl_yearly$time.series[, "remainder"] 
  
  # Save the stationary component in the matrix
  stationary_matrix[, i] <- deseasonalised
  
  
}


stationary_data <- as.data.frame(stationary_matrix)
colnames(stationary_data) <- colnames(POR_data)
head(stationary_data)

plot(stationary_data[,1], type = "l")
acf(stationary_data[,1], lag=200)

# Save the deseasonalised data
output_file <- "stationary_data.csv"
write.csv(stationary_data, file = output_file, row.names = FALSE)


# Plot the acf boxplots and the time series in Lisbon
lag_length <- 48

acf_list <- lapply(1:24, function(i) stats::acf(stationary_data[, i], lag.max = lag_length, plot = FALSE))

acf_data <- data.frame(
  Lag = rep(0:lag_length, times = 24),  
  ACF_Value = unlist(lapply(acf_list, function(acf_result) acf_result$acf)),
  Column = rep(1:24, each = (lag_length+1)) 
)

ggplot(acf_data, aes(x = factor(Lag), y = ACF_Value)) +
  geom_boxplot() +
  labs(title = "", x = "Lag", y = "ACF") +
  scale_x_discrete(breaks = seq(0, lag_length, by = 24))+ 
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-conf_limit, conf_limit), 
             linetype = "dotted", color = "blue")+
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  ) 

ggsave("empacfs.eps", width = 20, height = 20, units = "cm")

# Plot Lisbon time series

Time <- data[, 1]
Time <- as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S")

Lisbon_data <- stationary_data[, 22]

time_series_data <- data.frame(Time, Lisbon_data)

ggplot(time_series_data, aes(x = Time, y = Lisbon_data)) +
  geom_line(color = "black", size = 1) +
  labs( x = "Time", y = "Lisbon time series data (des.)") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
ggsave("Lisbondata.eps", width = 20, height = 20, units = "cm")

##################Heatmap of stationary data
POR_data_combined <- cbind(Time = Time_ad, stationary_data)

POR_data_long <- melt(POR_data_combined, id.vars = "Time")

POR_data_long$variable <- as.numeric(gsub("X", "", POR_data_long$variable))

breaks <- seq(1, max(POR_data_long$variable), by = 4)


ggplot(POR_data_long, aes(x = variable, y = Time, fill = value)) +
  geom_tile() +
  scale_x_continuous(breaks = breaks) +  
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),    
    values = scales::rescale(c(-1, 0, 1)), 
    limits = c(-1, 1),                     
    name = "Value"
  ) +
  labs(
    x = "Component",
    y = "Time"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
ggsave("heatmap.eps", width = 20, height = 20, units = "cm")


######################
# Fit the Grap supOU model
# Read in data
my_data <- read.csv(file = "stationary_data.csv")

#Read in adjacency matrix
Port_adjacency <- apply(as.matrix(read.csv(file = "adj_mat.csv")),2, as.numeric)

#Column normalisation of adjacency matrix
A_norm <- col_normalise_adjacency_matrix(Port_adjacency)
A_norm

X<-as.matrix(my_data)

# Compute the parameters for different lags
lag_vector <- seq(5, 100, by=1)
l <- length(lag_vector)
alpha_vector <- numeric(l)
c_vector <- numeric(l)

ev_all <- find_max_eigen(X, 100)
max_real_eigenvalues_all <- ev_all$max_real_eigenvalues

for (i in 1:l) {
  print(i)
  max_lag <- lag_vector[i]
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
    lower = c(-0.9999, 1),  # Constraints: c > -1, alpha > 1
    upper = c(0.9999, 100)  # Constraints: c < 1, alpha > 1
  )
  print(fit$par)
  c_vector[i] <- fit$par[1]
  alpha_vector[i] <- fit$par[2]
}

# Prepare data for ggplot
data <- data.frame(Lag = lag_vector, alpha = alpha_vector, c = c_vector)
data_long <- pivot_longer(data, cols = c("alpha", "c"), names_to = "Variable", values_to = "Value")

ggplot(data_long, aes(x = Lag, y = Value, color = Variable)) +
  geom_line(size = 1) +  
  labs(
    x = "Lag",
    y = "Values",
    color = "" 
  ) +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c(expression(alpha), "c")         
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
ggsave("alphac-vec.eps", width = 20, height = 20, units = "cm")


# Fit model for 40 lags
max_lag <- 40
ev <- find_max_eigen(X, max_lag)
max_real_eigenvalues <- ev$max_real_eigenvalues

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
  lower = c(-0.9999, 1),  # Constraints: c > -1, alpha > 1
  upper = c(0.9999, 100)  # Constraints: c < 1, alpha > 1
)
fit$par


# Compute function values
c <- fit$par[1] 
alpha <- fit$par[2]
y_values <- f(h_values, c, alpha, Delta)


# Plot the function
plot(
  h_values, y_values,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "h (Lag)",
  ylab = expression((1 + (1+c) * h * Delta)^(1 - alpha)),
  main = "Plot of the fitted eigenvalue function ",
  ylim = c(0, max(c(y_values, max_real_eigenvalues)))
)

grid()

points(
  h_values, max_real_eigenvalues,
  type = "o",
  col = "red",
  pch = 17
)

# Plot diagnostic plots for all eigenvalues


A_norm_eigen <- eigen(A_norm)$value



### OU case
initial_guess_OU <- 1  
fit_OU <- optim(
  par = initial_guess_OU,
  fn = objective_function_OU,
  h = h_values,
  Delta=Delta,
  observed = max_real_eigenvalues,
  method = "L-BFGS-B",
  lower = 0,  
  upper = Inf 
)
fit_OU$par


lambda_OU <- fit_OU$par

y_OU_values <- f_exp(h_values, lambda_OU, Delta)

#####
### 2OU case
initial_guess_2OU <- c(1,1,0.3)  
fit_2OU <- optim(
  par = initial_guess_2OU,
  fn = objective_function_2OU,
  h = h_values,
  Delta=Delta,
  observed = max_real_eigenvalues,
  method = "L-BFGS-B",
  lower = c(0,0,0),  
  upper = c(Inf,Inf,0.4999) 
)
fit_2OU$par


# Compute function values
lambda_2OU <- fit_2OU$par

y_2OU_values <- f_2exp(h_values, fit_2OU$par[1], fit_2OU$par[2], fit_2OU$par[3], Delta)


####Using ggplot
# Create a data frame for ggplot
plot_data <- data.frame(
  h = rep(h_values, 4),
  value = c(y_values, max_real_eigenvalues, y_OU_values, y_2OU_values),
  category = factor(
    rep(c("supOU - Gamma", "Empirical", "OU", "supOU - 2 Exp"), each = length(h_values)),
    levels = c("Empirical", "OU", "supOU - 2 Exp", "supOU - Gamma")
  )
)

colors <- c("Empirical" = "red", "OU" = "green", "supOU - 2 Exp" = "orange", "supOU - Gamma" = "blue")
shapes <- c("Empirical" = 16, "OU" = 17, "supOU - 2 Exp" = 0, "supOU - Gamma" = 5)

p <- ggplot(plot_data, aes(x = h, y = value, color = category, shape = category)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  labs(
    x = "h (Lag)",
    y = "Eigenvalue function",
    color = NULL,
    shape = NULL
  ) +
  theme(
    legend.position = "bottom",               
    legend.box = "horizontal",                
    legend.title = element_blank(),           
    legend.margin = margin(t = 5, b = 5, l = 5, r = 5), 
    legend.spacing.x = unit(0.3, "cm")        
  ) +
  guides(color = guide_legend(nrow = 2)) +    
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    text = element_text(size = 30),                 
    axis.title = element_text(size = 30),          
    axis.text = element_text(size = 30),           
    legend.text = element_text(size = 30),         
    legend.title = element_blank(),                
    legend.position = "bottom",                    
    legend.box = "horizontal",                     
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.caption = element_text(size = 30)         
  ) 

print(p)

ggsave("Fit.eps", plot = p, width = 20, height = 20, units = "cm")

#####Repeat for lag =100
max_lag <- 100
ev <- find_max_eigen(X, max_lag)
max_real_eigenvalues <- ev$max_real_eigenvalues

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
  upper = c(0.9999, 100)  
)
fit$par


# Compute function values
c <- fit$par[1] 
alpha <- fit$par[2]
y_values <- f(h_values, c, alpha, Delta)


# OU case
initial_guess_OU <- 1  
fit_OU <- optim(
  par = initial_guess_OU,
  fn = objective_function_OU,
  h = h_values,
  Delta=Delta,
  observed = max_real_eigenvalues,
  method = "L-BFGS-B",
  lower = 0,  
  upper = Inf 
)
fit_OU$par


lambda_OU <- fit_OU$par

y_OU_values <- f_exp(h_values, lambda_OU, Delta)

#####
# 2OU case
initial_guess_2OU <- c(1,1,0.3)  
fit_2OU <- optim(
  par = initial_guess_2OU,
  fn = objective_function_2OU,
  h = h_values,
  Delta=Delta,
  observed = max_real_eigenvalues,
  method = "L-BFGS-B",
  lower = c(0,0,0),  
  upper = c(Inf,Inf,0.4999) 
)
fit_2OU$par


lambda_2OU <- fit_2OU$par

y_2OU_values <- f_2exp(h_values, fit_2OU$par[1], fit_2OU$par[2], fit_2OU$par[3], Delta)




# Create a data frame for ggplot
plot_data <- data.frame(
  h = rep(h_values, 4),
  value = c(y_values, max_real_eigenvalues, y_OU_values, y_2OU_values),
  category = factor(
    rep(c("supOU - Gamma", "Empirical", "OU", "supOU - 2 Exp"), each = length(h_values)),
    levels = c("Empirical", "OU", "supOU - 2 Exp", "supOU - Gamma")
  )
)

colors <- c("Empirical" = "red", "OU" = "green", "supOU - 2 Exp" = "orange", "supOU - Gamma" = "blue")
shapes <- c("Empirical" = 16, "OU" = 17, "supOU - 2 Exp" = 0, "supOU - Gamma" = 5)

p <- ggplot(plot_data, aes(x = h, y = value, color = category, shape = category)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  labs(
    x = "h (Lag)",
    y = "Eigenvalue function",
    color = NULL,
    shape = NULL
  ) +
  theme(
    legend.position = "bottom",               
    legend.box = "horizontal",                
    legend.title = element_blank(),           
    legend.margin = margin(t = 5, b = 5, l = 5, r = 5), 
    legend.spacing.x = unit(0.3, "cm")        
  ) +
  guides(color = guide_legend(nrow = 2)) +    
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    text = element_text(size = 30),           
    axis.title = element_text(size = 30),     
    axis.text = element_text(size = 30),      
    legend.text = element_text(size = 30),    
    legend.title = element_blank(),         
    legend.position = "bottom",             
    legend.box = "horizontal",              
    plot.title = element_text(size = 30, hjust = 0.5), 
    plot.caption = element_text(size = 30)         
  ) 

print(p)

ggsave("Fit100.eps", plot = p, width = 20, height = 20, units = "cm")





###################################
# Recover the mean and variance matrix

max_lag <- 40
ev <- find_max_eigen(X, max_lag)
max_real_eigenvalues <- ev$max_real_eigenvalues

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
  upper = c(0.9999, 100)  
)
fit$par


c <- fit$par[1] 
alpha <- fit$par[2]

K <- create_K(c, A_norm)

mean <- (1-alpha)*K%*%as.numeric(colMeans(X))
var <- (1-alpha)*(K%*%cov(X)+cov(X)%*%t(K))


# Plot the heatmap for the mean vector
mean_data <- data.frame(Dimension = 1:24, Mean = mean)

p <- ggplot(mean_data, aes(x = Dimension, y = Mean, fill = Mean)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.0004, 0.031), name = "Mean") +
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

ggsave("estmean.eps", plot = p, width = 20, height = 20, units = "cm")

# Plot the heatmap for the variance-covariance matrix
var_data <- melt(var)
colnames(var_data) <- c("Row", "Column", "Value")
var_data$Row <- factor(var_data$Row)
var_data$Column <- factor(var_data$Column)

min(var_data$Value)
max(var_data$Value)

p<-ggplot(var_data, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.0004, 0.031), name = "Cov") +
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
  scale_y_discrete(breaks = paste0("X", seq(4, 24, by = 4)), labels = seq(4, 24, by = 4))

  
p
ggsave("estvar.eps",  plot=p, width = 20, height = 20, units = "cm")

# Add data covariance as comparison:

var_data <- melt(cov(X))
colnames(var_data) <- c("Row", "Column", "Value")

var_data$Row <- factor(var_data$Row)
var_data$Column <- factor(var_data$Column)

min(var_data$Value)
max(var_data$Value)

p<-ggplot(var_data, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.0004, 0.031), name = "Cov") +
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
  scale_x_discrete(breaks = paste0("X", seq(4, 24, by = 4)), labels = as.character(seq(4, 24, by = 4))) +
  scale_y_discrete(breaks = paste0("X", seq(4, 24, by = 4)), labels = seq(4, 24, by = 4))


p
ggsave("empvar.eps",  plot=p, width = 20, height = 20, units = "cm")


# Plot normalised A

var_data <- melt(A_norm)
colnames(var_data) <- c("Row", "Column", "Value")
var_data$Row <- factor(var_data$Row)
var_data$Column <- factor(var_data$Column)

min(var_data$Value)
max(var_data$Value)

p<-ggplot(var_data, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "black",  limits = c(0,1), name = expression(bar(A))) +
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
  scale_y_discrete(breaks = seq(4, 24, by = 4), labels = seq(4, 24, by = 4))


p
ggsave("normA.eps",  plot=p, width = 20, height = 20, units = "cm")

