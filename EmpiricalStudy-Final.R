# ============================================================
# EmpiricalStudy-Final.R
#
# Empirical application of the graph supOU model to wind
# capacity factors (WCFs) at 24 nodes in Portugal, using
# data from the RE-Europe dataset (Jensen & Pinson, 2017).
#
# This script reproduces the empirical study in Section 5 of:
#
#   Mehta, S. and Veraart, A. E. D. (2026).
#   "Statistical inference for Lévy-driven graph supOU processes:
#    From short- to long-memory in high-dimensional time series."
#   arXiv:2502.08838 [stat.ME] (preprint).
#   https://arxiv.org/abs/2502.08838
#
# Required input files (must be in the working directory):
#   wind_signal_COSMO.csv  -- hourly WCF time series, N=26304 obs,
#                             d=24 nodes; column 1 = timestamps,
#                             columns 2-25 = node WCFs
#   adj_mat.csv            -- 24x24 adjacency matrix for the
#                             Portuguese subnetwork (Figure 4a)
#
# Required R files (must be sourced before running):
#   GraphSupOU-Functions-Final.R  -- defines col_normalise_adjacency_matrix,
#                                    find_max_eigen, objective_function,
#                                    objective_function_OU, objective_function_2OU,
#                                    f, f_exp, f_2exp, create_K
#
# Output files (EPS):
#   heatmap-original.eps   -- heatmap of raw WCF time series (Figure 5a)
#   empacfs-original.eps   -- ACF boxplots of raw data         (Figure 5b)
#   Lisbondata-original.eps-- Lisbon node time series (raw)    (Figure 5c)
#   heatmap.eps            -- heatmap of deseasonalised data   (Figure 5d)
#   empacfs.eps            -- ACF boxplots after deseasonalising(Figure 5e)
#   Lisbondata.eps         -- Lisbon node time series (des.)   (Figure 5f)
#   alphac-vec.eps         -- alpha-hat and c-hat vs N*        (Figure 6a)
#   Fit.eps                -- fitted eigenvalue curves, N*=40  (Figure 6b)
#   Fit100.eps             -- fitted eigenvalue curves, N*=100 (Figure 6c)
#   estmean.eps            -- estimated Levy basis mean        (Figure 7a)
#   estvar.eps             -- estimated Levy basis variance     (Figure 7b)
#   empvar.eps             -- empirical covariance Var(X)      (Figure 7c)
#   normA.eps              -- column-normalised adjacency matrix(Figure 4b)
#   stationary_data.csv    -- deseasonalised and detrended WCFs
#
# Dependencies: forecast, ggplot2, reshape2, tidyr, ggrastr
# Minimum requirements: R >= 4.1.0
# Note: GraphSupOU-Functions-Final.R must be sourced first.
# ============================================================
library(forecast)
library(ggplot2)
library(reshape2)
library(tidyr)
library(ggrastr)

# ── Load raw data ─────────────────────────────────────────────────────────────
# wind_signal_COSMO.csv: N=26304 hourly observations (2012-2014),
# d=24 nodes in Portugal. Column 1 = timestamp, columns 2-25 = WCFs.
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


# ── Exploratory plots: raw data ───────────────────────────────────────────────
# Figures 5a, 5b, 5c in the paper.

# Prepare the data for plotting
Time_ad <- as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S")

POR_data_combined <- cbind(Time = Time_ad, POR_data)

POR_data_long <- melt(POR_data_combined, id.vars = "Time")

POR_data_long$variable <- as.numeric(gsub("X", "", POR_data_long$variable))

breaks <- seq(1, max(POR_data_long$variable), by = 4)

# Figure 5a: heatmap of raw 24-dimensional WCF time series
ggplot(POR_data_long, aes(x = variable, y = Time, fill = value)) +
  #geom_tile() +
  rasterise(geom_tile(), dpi = 150) +   # rasterize only the tiles
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
ggsave("heatmap-original.eps", width = 20, height = 20, units = "cm", dpi=150)


# ── Figure 5b: ACF boxplots of raw data ──────────────────────────────────────
# Boxplots of ACFs across all 24 nodes at lags 0, ..., 48.
# The dotted blue lines show the 95% confidence band +/- 1.96/sqrt(N).
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

# ── Figure 5c: Lisbon (node 22) raw time series ───────────────────────────────
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


# ── Deseasonalisation and detrending ─────────────────────────────────────────
# Raw WCFs exhibit daily (period = 24h) and yearly (period = 24*365h) seasonality
# and a mild trend. We remove both using two-pass STL decomposition:
#   Pass 1: remove daily seasonality and retain (remainder + trend)
#   Pass 2: remove yearly seasonality from the Pass 1 output
# The final "remainder" component is stored as the stationary series.
# See Section 5 of the paper for discussion.

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

# Save the deseasonalised data for downstream use
output_file <- "stationary_data.csv"
write.csv(stationary_data, file = output_file, row.names = FALSE)


# ── Exploratory plots: deseasonalised data ────────────────────────────────────
# Figures 5d, 5e, 5f in the paper.

# Figure 5e: ACF boxplots of deseasonalised data
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

# Figure 5f: Lisbon (node 22) deseasonalised time series
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

# Figure 5d: heatmap of deseasonalised data
POR_data_combined <- cbind(Time = Time_ad, stationary_data)

POR_data_long <- melt(POR_data_combined, id.vars = "Time")

POR_data_long$variable <- as.numeric(gsub("X", "", POR_data_long$variable))

breaks <- seq(1, max(POR_data_long$variable), by = 4)


ggplot(POR_data_long, aes(x = variable, y = Time, fill = value)) +
  #geom_tile() +
  rasterise(geom_tile(), dpi = 150) +   # rasterize only the tiles
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
ggsave("heatmap.eps", width = 20, height = 20, units = "cm", dpi=150)


# ── Fit the graph supOU model ─────────────────────────────────────────────────
# Estimation follows the two-step procedure of Section 3.1:
#   Step 1: minimise L(alpha, c) over the empirical max real eigenvalues
#           of R_hat(h) = cov_hat(h) %*% cov_hat(0)^{-1} (equation (14))
#   Step 2: recover mu_L and sigma^2_L from sample mean and covariance

# Read in data
my_data <- read.csv(file = "stationary_data.csv")

#Read in adjacency matrix
Port_adjacency <- apply(as.matrix(read.csv(file = "adj_mat.csv")),2, as.numeric)

# Column normalisation of adjacency matrix (see Section 2.3 and Figure 4b)
A_norm <- col_normalise_adjacency_matrix(Port_adjacency)
A_norm

X<-as.matrix(my_data)

# ── Figure 6a: alpha-hat and c-hat as a function of N* ───────────────────────
# Sweep N* from 5 to 100 lags. For each N*, re-estimate (alpha, c) by
# minimising the Gamma supOU loss function (14) via L-BFGS-B with
# box constraints: c in (-1, 1), alpha > 1.
# The estimates stabilise around N* = 40; see Section 5.

# Compute the parameters for different lags
lag_vector <- seq(5, 100, by=1)
l <- length(lag_vector)
alpha_vector <- numeric(l)
c_vector <- numeric(l)

# Precompute eigenvalues for all lags up to 100 once (reused in the loop)
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
    x = bquote(N^"*"),
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


# ── Figure 6b: fitted eigenvalue curves at N* = 40 ───────────────────────────
# Three models are compared against the empirical max real eigenvalue sequence:
#   - Graph supOU with Gamma(alpha, 1) kernel  [equation (12)]
#   - Graph OU (single exponential decay)
#   - Graph supOU with sum of two exponentials [equation (11)]

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


# Quick diagnostic plot (base R) comparing fitted vs empirical eigenvalue curve
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

# ── Figure 6c: fitted eigenvalue curves at N* = 100 ──────────────────────────
# Repeat the three-model comparison using N* = 100 lags to confirm
# robustness of the estimates to the choice of N*; see Section 5.

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


# ── Step 2: recover Lévy basis mean and variance ──────────────────────────────
# Using the N*=40 estimates of (alpha, c), recover the Lévy basis parameters:
#   mu_L  = (1 - alpha) * K(c) %*% X_bar
#   sigma^2_L = (1 - alpha) * (K(c) %*% Cov(X) + Cov(X) %*% K(c)^T)
# See equation (9) and Proposition 3.2 in the paper.
# Figures 7a, 7b, 7c in the paper.

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

# Construct drift matrix K(c) = -I - c * A_bar^T (equation (5))
K <- create_K(c, A_norm)

# Estimate Levy basis mean and variance (equation (9), Proposition 3.2)
mean <- (1-alpha)*K%*%as.numeric(colMeans(X))
var <- (1-alpha)*(K%*%cov(X)+cov(X)%*%t(K))


# ── Figure 7a: estimated Lévy basis mean ─────────────────────────────────────
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

# ── Figure 7b: estimated Lévy basis variance matrix ──────────────────────────
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

# ── Figure 7c: empirical covariance Var(X) for comparison ────────────────────
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


# ── Figure 4b: column-normalised adjacency matrix ─────────────────────────────

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