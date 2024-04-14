##############
#
# Chapter 1
#
##############

library(ggplot2)
library(latex2exp)
library(tinytex)
library(ggthemes)
library(ggExtra)
library(reshape2)

#
#Brownian Motion
#
mo.bro.w.dist <- function(T = 1, # Time in years
                          n = 10, # Number of steps
                          N = 5){ # Number of sample paths
  dt <- T / n # Step size
  t.seq <- seq(0, T, by = dt) # Sequence of time points from zero to T in dt sized increments
  B <- matrix(NA, nrow = (n+1), ncol = N) # An empty (n+1) by N matrix
  for (j in 1:N) {
    xx <- c(0, rnorm(n, mean = 0, sqrt(dt))) # Generate n random normal values with variance T / n = dt and initial value zero
    yy <- cumsum(xx) # Brownian motion is the cumulative sum of generated values
    B[,j] <- yy # Update column j of matrix B with values of yy 
  }
  
# ggMarginal(), the function used to plot the graph requires a data frame input therefore we must convert the matrix B into a suitable data frame 
  Bt <- t(B) # Transpose matrix B
  rownames(Bt) <- paste("trial", seq(N), sep = "") # Set row names to to trial number ( Trial 1 to Trial N)
  colnames(Bt) <- paste("time", t.seq, sep = "") # Set column headers to time sequence values
  
  dat <- as.data.frame(Bt) # Convert matrix Bt to a data frame 
  dat$trial <- rownames(dat) # Add additional column labeled trial to the end of the data frame
  mdat <- melt(dat, id.vars="trial") # Reshape and elongate the data frame to time order with each time block ordered by trial
  mdat$time <- as.numeric(gsub("time", "", mdat$variable)) # Add additional column labeled time with just numeric values for time
  
  
  p <- ggplot(mdat, aes(x = time, y = value, group = trial)) + # Plot value (cumulative sum of random values) against time from the new data frame mdat
    theme_bw() + # Black and white background
    theme(panel.grid=element_blank()) + # Remove grid lines
    geom_line(linewidth = 0.2, alpha = 0.5, aes(color = trial))+ # Join data points with a line and colour each trial independently
    theme(legend.position = "none",  # No legend is needed as each line is an independent path of Brownian motion and knowing which trial is which is not necassary
          plot.caption = element_text(hjust = 0.5, size = 15), # Set caption size and position
          axis.title.x = element_text(size = 17), # Set x axis title font size
          axis.title.y = element_text(size = 17), # set y axis title font size
          axis.text.x = element_text(size = 14), # Set x axis values font size
          axis.text.y = element_text(size = 14))+ # Set y axis values font size
    geom_point(size = 0.0000000001, aes(color = trial))+ # Set point size of data to very small to remove noise when a large number of trials are being used
    geom_hline(yintercept = 0) +  # Horizontal line at y=0
    labs(x = TeX("Time $t$ (Years)"),  # x axis label
         y = TeX("$B_t$ ~ $N(0,dt)$")) # y axis label
  ggMarginal(p, type = "density", margins = 'y') # plot p along with a density function along the y axis
}
mo.bro.w.dist(n = 500, T = 1, sigma = 1, N = 50)

#
#Brownian Motion with Drift
#
mo.bro.drift.w.dist <- function(T = 1, # Time in years
                                sigma = 1, # Volatility
                                mu = 12, # Drift term
                                n = 500, # Number of steps  
                                N = 50, # Number of sample paths
                                do_plot = FALSE){ 
  dt <- T / n # Step size
  t.seq <- seq(0,T,by = dt) # Sequence of time points from zero to T in dt sized increments
  B_d <- matrix(NA, nrow = (n+1), ncol = N) # An empty (n+1) by N matrix
  for (j in 1:N) {
    xx <- c(0, rnorm(n, mean = 0, sqrt(dt))) # Generate n random normal values with variance T / n = dt and initial value zero
    yy <- cumsum(xx) # Brownian motion is the cumulative sum of generated values
    Bd <- mu * t.seq + sigma * yy # Brownian motion with drift value = mu * t + sigma * B_t
    B_d[,j] <- Bd # Update column j of matrix B with values of Bd 
  }
  
# Data frame and graph are formed as in mo.bro.w.dist (line 27 to 51)
  Bt <- t(B_d)
  rownames(Bt) <- paste("trial", seq(N), sep = "")
  colnames(Bt) <- paste("time", t.seq, sep = "")
  
  dat <- as.data.frame(Bt)
  dat$trial <- rownames(dat)
  mdat <- melt(dat, id.vars = "trial")
  mdat$time <- as.numeric(gsub("time", "", mdat$variable))
  
  leg.txt <- TeX(sprintf("$X_t = mu t =$ %s $t$", signif(mu, 2)))
  
  p <- ggplot(mdat, aes(x = time, y = value, group = trial)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          plot.caption = element_text(hjust = 0.5, size = 15),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    geom_line(linewidth = 0.2, alpha = 0.5, aes(color = trial)) +
    theme(legend.position = "none") +
    geom_point(size = 0.0000000001, aes(color = trial))+
    # Add the mean drift term Xt = mu*t 
    geom_abline(slope = mu, 
                intercept = 0,
                color="red") + 
    labs(x = TeX("Time $t$ (Years)"), 
         y = TeX("$X_t = mu t + sigma B_t$"),
         caption = TeX(sprintf("$sigma$ = %s, $mu$ = %s", sigma, mu)))
  
  
  if(do_plot){
    ggMarginal(p, type = "density", margins = 'y') # If do_plot TRUE, plot density function
  }
  else{
    p # Else just plot p
  }  
}

#Positive drift
mo.bro.drift.w.dist(n = 500, T = 1, sigma = 1, mu = 12, N = 50)

# Negative drift
mo.bro.drift.w.dist(n = 500, T = 1, sigma = 1, mu = -12, N = 50)

#
# Geometric Brownian motion
#
geo.mo.brow.dist <- function(T = 1, # Time in years
                             sigma = 1, # Volatility
                             mu = 0.5, # Drift term
                             X_0 = 24, # Initial value
                             n = 500, # Number of steps
                             N = 50, # Number of sample paths
                             do_plot = TRUE){
  dt <- T / n # Step size                                   
  t.seq <- seq(0,T,by = dt) # Sequence of time points from zero to T in dt sized increments
  X <- matrix(NA, nrow = (n+1), ncol = N)  # An empty (n+1) by N matrix
  for (j in 1:N){
    xx <- c(0, rnorm(n, mean = 0, sqrt(dt))) # Generate n normal random variates with variance T / n = dt
    yy <- cumsum(xx) # Brownian Motion is the cumulative sum of the random variates
    zz <- X_0 * exp(mu * t.seq + sigma * yy) # Geometric Brownian motion 
    X[,j] <- zz  # Update column j of matrix X with values of zz
  }
  
# Data frame and graph are formed as in mo.bro.w.dist (line 20 to 44)
  Bt <- t(X)
  rownames(Bt) <- paste("trial", seq(N), sep = "")
  colnames(Bt) <- paste("time", t.seq, sep = "")
  
  dat <- as.data.frame(Bt)
  dat$trial <- rownames(dat)
  mdat <- melt(dat, id.vars = "trial")
  mdat$time <- as.numeric(gsub("time", "", mdat$variable))
  
  
  p <- ggplot(mdat, aes(x = time, y = value, group = trial)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          plot.caption = element_text(hjust = 0.5, size = 15),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    geom_line(linewidth = 0.2, alpha = 0.5, aes(color = trial)) +
    theme(legend.position = "none") +
    geom_point(size = 0.0000000001, aes(color = trial)) + 
    labs(x = TeX("Time $t$ (Years)"), 
         y = TeX("$X_t = X_0 e^{mu t + sigma B_t}$"),
         caption = TeX(sprintf("$sigma$ = %s, $mu$ = %s, $X_0$ = %s", signif(sigma,4), signif(mu,4), signif(X_0,4))))
  
  if(do_plot){
    ggMarginal(p, type = "density", margins = 'y') # If do_plot TRUE, plot density function
  }
  else{
    p # Else just plot p
  }
}

geo.mo.brow.dist(N = 50, n = 500, T = 1, sigma = 1, mu = 0.5, X_0 = 24)



##############
#
# Chapter 2
#
##############

library(car)
library(quantmod)
library(MASS)

# We want our step size to be daily so for a 252 trading days per year we have the delta value 1/252.
delta <- 1 / 252

#
# BMW
#

bmw_1 <- getSymbols(Symbols = "BMW.DE", # Get data for BMW
                    src = "yahoo", # Data provider 
                    auto.assign = FALSE, # Allows data to be assigned directly to the bmw_1 object and only works for downloading data on a single stock
                    return.class = "xts") # Return data as a time series
bmw <- na.omit(bmw_1)  # This removes any NA values from the time series.

# Extract the closing prices
BMW_Close <- bmw[760:3296,"BMW.DE.Close"]

# Number of data points, excluding first
m_1 <- length(BMW_Close) - 1

# Extract the values themselves
BMW_Close_values <- drop(coredata(BMW_Close))

#
plot(BMW_Close_values,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = "BMW Closing Values(€)")

# Difference of logs
BMW_diff_log <- diff(log(BMW_Close_values))


# Estimate the parameters of geometric Brownian motion
#
# Estimate mu
mu_hat_bmw <- mean(BMW_diff_log) / delta
#
# Estimate sigma^2
sigma_2_hat_bmw <- var(BMW_diff_log) / delta
#
# Estimate sigma
sigma_hat_bmw <- sqrt(sigma_2_hat_bmw)
#
c(mu_hat_bmw, sigma_2_hat_bmw, sigma_hat_bmw) # Returns values for mu, sigma and sigma^2
#
# Generate 5 sample paths of geometric Brownian motion for BMW share price data
geo.mo.brow.dist(T = m_1/252, # m_1 is daily data so needs dividing by 252 to get time in years
                 sigma = sigma_hat_bmw, 
                 mu = mu_hat_bmw, 
                 X_0 = BMW_Close_values[1],
                 n = m_1,
                 N = 10,
                 do_plot = FALSE)


#
# JP Morgan Chase & Co.
#

jpm <- getSymbols(Symbols = "JPM", # Get data for BMW
                  src = "yahoo", 
                  auto.assign = FALSE,
                  return.class = "xts")

# Extract the closing prices from 31.12.2009 to 31.12.2019
JPM_Close <- jpm[756:3272,"JPM.Close"]

# Number of data points, excluding first
m_2 <- length(JPM_Close) - 1

# Extract the values themselves
JPM_Close_values <- drop(coredata(JPM_Close))

# Difference of log share price values
JPM_diff_log <- diff(log(JPM_Close_values))

# Estimate the parameters of geometric Brownian motion
#
# Estimate mu
mu_hat_jpm <- mean(JPM_diff_log) / delta
#
# Estimate sigma^2
sigma_2_hat_jpm <- var(JPM_diff_log) / delta
#
# Estimate sigma
sigma_hat_jpm <- sqrt(sigma_2_hat_jpm)
#
c(mu_hat_jpm, sigma_2_hat_jpm, sigma_hat_jpm)
#

#
# Atlas Copco
#

A_C <- getSymbols(Symbols = "ATCO-A.ST", # Get data for BMW
                  src = "yahoo", 
                  auto.assign = FALSE,
                  return.class = "xts")

# Extract the closing prices from 31.12.2009 to 31.12.2019
A_C_Close <- A_C[753:3264,"ATCO-A.ST.Close"]

# Number of data points, excluding first
m_3 <- length(A_C_Close) - 1

# Extract the values themselves
A_C_Close_values <- drop(coredata(A_C_Close))

# Difference of logs
A_C_diff_log <- diff(log(A_C_Close_values))

# Estimate the parameters of geometric Brownian motion
#
# Estimate mu
mu_hat_A_C <- mean(A_C_diff_log) / delta
#
# Estimate sigma^2
sigma_2_hat_A_C <- var(A_C_diff_log) / delta
#
# Estimate sigma
sigma_hat_A_C <- sqrt(sigma_2_hat_A_C)
#
c(mu_hat_A_C, sigma_2_hat_A_C, sigma_hat_A_C)
#


#
# Difference of log share price
#


#Simulating 2520 random normal values
N <- 2520 # Equivalent to ten years of data
rand_norm <- rep(NA, N+1) 
N_seq <- seq(from = 0,
             to = N,
             length = N+1)
for (i in N_seq){
  rand_norm[i-1] <- rnorm(1, mean = 10, sd =0.05) 
}
# Take the difference of the log values
Norm_diff_log <- diff(log(rand_norm)) 


par(mfrow = c(2,2))
plot(BMW_diff_log,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = " ",
     main = "BMW",
     sub = "(a)")
#
plot(JPM_diff_log,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = " ",
     main = "JPMorgan Chase & Co.",
     sub = "(b)")
#
plot(A_C_diff_log,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = " ",
     main = "Atlas Copco",
     sub = "(c)")
#
plot(Norm_diff_log,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = " ",
     main = "Simulated Normal Values",
     sub = "(d)")

# Estimate the parameters of t-distribution

# BMW
#
# Fit a t distribution to the BMW share price data
BMW_t_fit <- fitdistr(BMW_diff_log, "t")
#
# Location Estimate
l_hat_bmw <- BMW_t_fit$estimate[1] 
#
# Scale Estimate
s_hat_bmw <- BMW_t_fit$estimate[2]
#
# Degrees of Freedom Estimate
nu_hat_bmw <- BMW_t_fit$estimate[3]

# JPM
#
# Fit a t distribution to the JPM share price data
JPM_t_fit <- fitdistr(JPM_diff_log, "t")
#
# Location Estimate
l_hat_jpm <- JPM_t_fit$estimate[1]
#
# Scale Estimate
s_hat_jpm <- JPM_t_fit$estimate[2]
#
# Degrees of Freedom Estimate
nu_hat_jpm <- JPM_t_fit$estimate[3]

# Atlas Copco
#
# Fit a t distribution to the JPM share price data
A_C_t_fit <- fitdistr(A_C_diff_log, "t")
#
# Location Estimate
l_hat_A_C <- A_C_t_fit$estimate[1]
#
# Scale Estimate
s_hat_A_C <- A_C_t_fit$estimate[2]
#
# Degrees of Freedom Estimate
nu_hat_A_C <- A_C_t_fit$estimate[3]

#
# QQ Plots
#

par(mfrow = c(3, 2))
#
#BMW
qqPlot(BMW_diff_log,
       distribution = "norm",
       main = "Normal Q-Q Plot of the Log Difference of BMW Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
qqPlot(BMW_diff_log,
       distribution = "t",
       df = nu_hat_bmw,
       main = "t Q-Q Plot of the Log Difference of BMW Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
# JPM
qqPlot(JPM_diff_log,
       distribution = "norm",
       main = "Normal Q-Q Plot of the Log Difference of JP Morgan Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
qqPlot(JPM_diff_log,
       distribution = "t",
       df = nu_hat_jpm,
       main = "t Q-Q Plot of the Log Difference of JP Morgan Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
#Atlas Copco
qqPlot(A_C_diff_log,
       distribution = "norm",
       main = "Normal Q-Q Plot of the Log Difference of Atlas Copco Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
qqPlot(A_C_diff_log,
       distribution = "t",
       df = nu_hat_A_C,main = "t Q-Q Plot of the Log Difference of Atlas Copco Stock Prices",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")


# 
# Share price model based on the t distribution
#

t.estimation.dist <- function(T = 1,# Time in years
                              s = 1, # Scale
                              loc = 0.5, # Location
                              X_0 = 24, # Initial value
                              nu = 3, # Degrees of freedom
                              n = 500, # Number of steps
                              N = 50, # Number of sample paths
                              do_plot = TRUE){
  dt <- T / n # Step size                                     
  t.seq <- seq(0, T, by = dt) # Sequence of time points from zero to T in dt sized increments
  X <- matrix(NA, nrow = (n+1), ncol = N) # An empty (n+1) by N matrix
  for (j in 1:N){
    xx <- c(0, rt(n, df=nu))   # Generate n t random variates with nu degrees of freedom
    yy <- cumsum(xx)    # Cumulative sum of the random variates
    zz <- X_0 * exp(n*loc + s*yy) # The t based model Xt = X0 e ^{n*l + s(t1+...+tN)} 
    X[,j] <- zz   # Update column j of matrix X with values of zz
  }
  
  # Data frame and graph are formed as in mo.bro.w.dist (line 27 to 51)  
  Bt <- t(X)
  rownames(Bt) <- paste("trial", seq(N), sep = "")
  colnames(Bt) <- paste("time", t.seq, sep = "")
  
  dat <- as.data.frame(Bt)
  dat$trial <- rownames(dat)
  mdat <- melt(dat, id.vars = "trial")
  mdat$time <- as.numeric(gsub("time", "", mdat$variable))
  
  
  p <- ggplot(mdat, aes(x = time, y = value, group = trial)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(linewidth = 0.2, alpha = 0.5, aes(color = trial)) +
    theme(legend.position = "none") +
    geom_point(size = 0.0000000001, aes(color = trial)) + 
    labs(x = "Time (Years)", 
         y = TeX("$X_t = X_0 e^{nl+s(t_{\nu,1}+t_{\nu,2}+...+t_{\nu,n})}$"))
  if(do_plot){
    ggMarginal(p, type = "density", margins = 'y')
  }
  else{
    p
  }
}


t.estimation.dist(T = m_1/252, 
                  s = s_hat_bmw,
                  loc = l_hat_bmw, 
                  X_0 = BMW_Close_values[1], 
                  nu = nu_hat_bmw,
                  n = 500,
                  N = 10,
                  do_plot = FALSE)


# Appendix A

#
# JPM
#

# Plot of JP Morgan Chase & Co. closing values from 2010 to 2019
plot(JPM_Close_values,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = "JPMorgan Chase Closing Values ($)")

# 5 simulated paths of JP Morgan Chase & Co. share price using geometric Brownian motion
geo.mo.brow.dist(T = m_1/252, 
                 sigma = sigma_hat_jpm,
                 mu = mu_hat_jpm,
                 X_0 = JPM_Close_values[1],
                 n = m_1,
                 N = 5,
                 do_plot = FALSE)



#
# Atlas Copco
#

# Plot of Atlas Copco closing values from 2010 to 2019
plot(A_C_Close_values,
     type = "l",
     xlab = TeX("Time (Days)"),
     ylab = "Atlas Copco Closing Values (kr)")

# 5 simulated paths of Atlas Copco share price using geometric Brownian motion
geo.mo.brow.dist(T = m_1/252,
                 sigma = sigma_hat_A_C,
                 mu = mu_hat_A_C, 
                 X_0 = A_C_Close_values[1],
                 n = m_1,
                 N = 5,
                 do_plot = FALSE)
#




##############
#
# Chapter 4
#
##############

# Monte Carlo Simulation for option price excluding hedging factor

get_price_Monte_Carlo_N <- function(T = 0.25, # Maturity in years
                                    sigma = 0.45, # Volatility
                                    mu = 1, # Drift term
                                    X_0 = 40, # Initial value
                                    K = 60, # Strike price 
                                    r = 0.06, # Risk free rate
                                    N = 1000){  # Number of Monte Carlo samples
  #
  B_T <- rnorm(N, 0, sqrt(T)) # We can do this because we know the distribution of B_T
  # This means that we don't have to simulate the complete sample path
  #
  X_T <- X_0 * exp(mu * T + sigma * B_T) # Share price value at time T
  #
  simulated_prices <- exp(-r * T) * pmax(X_T - K , 0) # Option Price discounted back to price at t=0
  #
  price_N <- mean(simulated_prices) # Average price of all simulated option prices
  #
  # Get the standard error
  #
  se_N <- sd(simulated_prices) / sqrt(N)
  #
  c(price_N, se_N)
}

get_price_Monte_Carlo_N(T = 1,
                        sigma = sigma_hat_bmw,
                        mu = mu_hat_bmw,
                        X_0 = BMW_Close_values[1],
                        K = BMW_Close_values[252], 
                        r = 0.005,  
                        N = 1000)


# Monte Carlo Simulation for option price including hedging factor
get_price_Monte_Carlo_H <- function(T = 0.25,
                                    sigma = 0.45,
                                    mu = 1, 
                                    X_0 = 60,
                                    K = 40, 
                                    r = 0.005, 
                                    N = 1000){
  #
  q <- (mu + 0.5 * sigma^2 - r) / sigma
  #
  # Get price
  #
  B_T <- rnorm(N, 0, sqrt(T)) # We can do this because we know the distribution of B_T
  # This means that we don't have to simulate the complete sample path
  #
  X_T <- X_0 * exp(mu * T + sigma * B_T) # Share price value at time T
  #
  simulated_prices_H <- exp(-r * T) * pmax(X_T - K , 0) * exp(-q * B_T - 0.5 * q^2 * T) # Option Price discounted back to price at t=0 including hedging factor
  #
  price_H <- mean(simulated_prices_H)  # Average price of all simulated option prices
  #
  # Get the standard error
  #
  se_H <- sd(simulated_prices_H) / sqrt(N)
  #
  c(price_H, se_H)
}

get_price_Monte_Carlo_H(T = 1,
                        sigma = sigma_hat_bmw,
                        mu = mu_hat_bmw,
                        X_0 = BMW_Close_values[1],
                        K = BMW_Close_values[252], 
                        r = 0.06,  
                        N = 1000)
#


BlackScholes <- function(T = 0.25, # Maturity in years
                         sigma = 0.45, # Volatility
                         X_0 = 20, # Initial value
                         K = 40, # Strike price
                         r = 0.005){ # Risk free rate
  #
  d1 <- (log(X_0 / K) + (r + sigma^2 / 2)*T) / (sigma * sqrt(T))
  #
  d2 <- d1 - sigma * sqrt(T)
  #
  value <- X_0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  #
  value
}

##############
#
# Chapter 5
#
##############
#
# Function to price option using normal increments
#
option_normal_increments <- function(T = 0.25, # Years 
                                     sigma = 0.45,
                                     mu = 0.001, 
                                     X_0 = 60,
                                     K = 40, 
                                     r = 0.06,
                                     delta = 1 / 252,
                                     N = 10000){
  #
  # Number of points at which to simulate process
  #
  n <- T / delta
  #
  sigma_2 <- sigma^2
  #
  # Sequence of n time values (after t = 0)
  #
  t_seq <- delta * (1:n)
  #
  # Estimate q
  #
  q <- (mu + 0.5 * sigma_2 - r) / sigma
  #
  # Save the naive prices, the correction term for hedging and the corrected price
  Price_naive_normal <- rep(NA, N)
  Correction_term_normal <- rep(NA, N)
  Price_corrected_normal <- rep(NA, N)
  #
  # Now generate N prices
  #
  for(j in 1:N){
    #
    # Generate Brownian motion
    #
    Z <- rnorm(n)
    B <- sqrt(delta) * cumsum(Z)
    #
    # Final value of B
    #
    B_T <- B[length(B)]
    #
    # Get the share prices (after t = 0)
    #
    X <- X_0 * exp(mu * t_seq + sigma * B)
    #
    # Final value
    #
    X_T <- X[length(X)]
    #
    # Naive price
    #
    Price_naive_normal[j] <- exp(-r * T) * max(X_T - K, 0)
    #
    # Correction term
    #
    Correction_term_normal[j] <- exp(-q * B_T  - 0.5 * q^2 * T)
    #
    # Corrected price
    #
    Price_corrected_normal[j] <- Price_naive_normal[j] * Correction_term_normal[j]
    #
  } # End of loop
  #
  option_price_normal <- mean(Price_corrected_normal)
  #
  return(option_price_normal)
  #
} # End of function
#

#
# Function to price option using t distribution
#
option_t_increments <- function(T = 0.25, # Maturity in years
                                sigma = 0.45,
                                mu = 0.001, 
                                X_0 = 43.8,
                                K = 45, 
                                r = 0.06,  
                                nu_hat = 3.44, # Degrees of freedom
                                delta = 1 / 252,
                                N = 10000,
                                do_trim = FALSE,
                                trim = 0.01){
  #
  # Number of points at which to simulate process
  #
  n <- T / delta
  #
  l <- mu * delta
  s <- sigma * sqrt(delta)
  #
  sigma_2 <- sigma^2
  #
  # Sequence of n time values (after t = 0)
  #
  t_seq <- delta * (1:n)
  #
  # Estimate q
  #
  q_hat <- (mu + 0.5 * sigma_2 - r) / sigma
  #
  # Save the naive prices, the correction term for hedging and the corrected price
  Price_naive_t <- rep(NA, N)
  Correction_term_t <- rep(NA, N)
  Price_corrected_t <- rep(NA, N)
  #
  # Now generate N prices
  #
  for(j in 1:N){
    #
    # Generate t random variates
    #
    T_r <- rt(n, nu_hat)
    #
    # Cumulative sum
    #
    S <- cumsum(T_r)
    #
    # Get the share prices (after t = 0)
    #
    X <- X_0 * exp(l * c(1:n) + s * S)
    #
    # Final value
    #
    X_T <- X[length(X)]
    #
    # Price (which is saved) [using same final term as in the normal increments case]
    #
    # We need B_T
    #
    # Transform T_r to standard normal
    #
    Z <- qnorm(pt(T_r, nu_hat))
    #
    B <- sqrt(delta) * cumsum(Z)
    #
    B_T <- B[length(B)]
    #
    # Naive price
    #
    Price_naive_t[j] <- exp(-r * T) * max(X_T - K, 0)
    #
    # Correction term
    #
    Correction_term_t[j] <- exp(-q_hat * B_T  - 0.5 * q_hat^2 * T)
    #
    # Corrected price
    #
    Price_corrected_t[j] <- Price_naive_t[j] * Correction_term_t[j]
    #
  }
  #
  if(do_trim){
    #
    option_price_t <- mean(Price_corrected_t, 
                           trim = trim) } else {
                             option_price_t <- mean(Price_corrected_t)
                             #
                           }
  #
  return(option_price_t)
  #
} # End of function



# Setting constant variables for option price
N_seq <- 20 # Number of data points
#
T <- 1 # Maturity in years
#
sigma <- sigma_hat_bmw # volatility of BMW share price
#
mu <- mu_hat_bmw # mean of BMW share price
#
X_0 <- BMW_Close_values[m_1-249] # Initial value of BMW share price on 02-01-2019
#
K <- BMW_Close_values[m_1+1] # Strike price (BMW share price value on 30-12-2019)
#
r <- 0.005 # Risk free rate in 2019, 0.5%
#
nu_hat <- nu_hat_bmw # Estimated degrees of freedom of t model fitted to BMW share price
#
delta <- 1/252
#
N <- 10000 # Number of Monte Carlo simulations



#
# How option price depends on interest Rate
#
# Sequence of 20 interest rate values from 0.1% to 5%
r_seq <- seq(from = 0.001,
             to = 0.05,
             length = N_seq)
#
#
price_t_nu_hat_seq <- rep(NA, N_seq)
#
price_t_3_seq <- rep(NA, N_seq)
#
price_t_10_seq <- rep(NA, N_seq)
#
price_t_500_seq <- rep(NA, N_seq)
#
price_t_normal_seq <- rep(NA, N_seq)
#
price_normal_seq <- rep(NA, N_seq)
#
Black_Scholes <- rep(NA, N_seq)

#
for(k in 1:N_seq){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T, 
                                               sigma = sigma,
                                               mu = mu,
                                               X_0 = X_0,
                                               K = K, 
                                               r = r_seq[k],
                                               nu_hat = nu_hat,
                                               delta = 1 / 252,
                                               N = N,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T,
                                          sigma = sigma,
                                          mu = mu,
                                          X_0 = X_0,
                                          K = K, 
                                          r = r_seq[k], 
                                          nu_hat = 3, 
                                          delta = delta,
                                          N = N,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T, 
                                           sigma = sigma,
                                           mu = mu,
                                           X_0 = X_0,
                                           K = K, 
                                           r = r_seq[k], 
                                           nu_hat = 10, 
                                           delta = delta,
                                           N = N,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T, 
                                            sigma = sigma,
                                            mu = mu,
                                            X_0 = X_0,
                                            K = K, 
                                            r = r_seq[k], 
                                            nu_hat = 500, 
                                            delta = delta,
                                            N = N,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T,
                                                  sigma = sigma,
                                                  mu = mu,
                                                  X_0 = X_0,
                                                  K = K,
                                                  r = r_seq[k],
                                                  delta = delta,
                                                  N = N)
  #
  Black_Scholes[k] <- BlackScholes(T = T,
                                   sigma = sigma,
                                   X_0 = X_0,
                                   K = K,
                                   r = r_seq[k])
}
#
#
par(mfrow = c(1,1))

#Legend
leg.txt.r <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(r_seq,  # interest rate on x axis
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b", # Points connected by a line
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Interest Rate"), # x axis label
        ylab = "Option Price", # y axis label
        main = TeX("The Dependence of Option Price on Interest Rate"), # Title
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $T$ = %s", signif(X_0,4), signif(K,4), signif(sigma,3), signif(T,1))), # Subtitle stating values of constant variables
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.r, pch=4, title= "Distribution", # Legend size and location
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)


#
# How option price depends on varying volatility
#
#######################

# Sequence of 20 volatility values from 1% to 40%
sigma_seq <- seq(from = 0.01,
                 to = 0.4,
                 length = N_seq)
#
price_t_nu_hat_seq <- rep(NA, N_seq)
#
price_t_3_seq <- rep(NA, N_seq)
#
price_t_10_seq <- rep(NA, N_seq)
#
price_t_500_seq <- rep(NA, N_seq)
#
price_normal_seq <- rep(NA, N_seq)
#
Black_Scholes <- rep(NA, N_seq)
#

#
for(k in 1:N_seq){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T, 
                                               sigma = sigma_seq[k],
                                               mu = mu,
                                               X_0 = X_0,
                                               K = K, 
                                               r = r, 
                                               nu_hat = nu_hat, 
                                               delta = delta,
                                               N = N,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T, 
                                          sigma = sigma_seq[k],
                                          mu = mu,
                                          X_0 = X_0,
                                          K = K, 
                                          r = r, 
                                          nu_hat = 3, # Degrees of freedom
                                          delta = delta,
                                          N = N,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T, 
                                           sigma = sigma_seq[k],
                                           mu = mu,
                                           X_0 = X_0,
                                           K = K, 
                                           r = r, 
                                           nu_hat = 10, # Degrees of freedom
                                           delta = delta,
                                           N = N,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T,
                                            sigma = sigma_seq[k],
                                            mu = mu, 
                                            X_0 = X_0,
                                            K = K, 
                                            r = r, 
                                            nu_hat = 500, # Degrees of freedom
                                            delta = delta,
                                            N = N,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T,
                                                  sigma = sigma_seq[k],
                                                  mu = mu,
                                                  X_0 = X_0,
                                                  K = K,
                                                  r = r,
                                                  delta = delta,
                                                  N = N)
  #
  Black_Scholes[k] <- BlackScholes(T = T,
                                   sigma = sigma_seq[k],
                                   X_0 = X_0,
                                   K = K,
                                   r = r)
}
#
#
par(mfrow = c(1,1))

#Legend
leg.txt.sigma <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes")
#
matplot(sigma_seq, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$\\sigma$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $\\sigma$"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $r$ = %s, $T$ = %s", signif(X_0,4), signif(K,4), signif(r,3), signif(T,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.sigma, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)


#
# Varying Strike Price
#
#######################

# Sequence of 20 possible strike prices from the initial value X_0 to X_0+10
K_seq <- seq(from = X_0,
             to = X_0+10,
             length = N_seq)
#
price_t_nu_hat_seq <- rep(NA, N_seq)
#
price_t_3_seq <- rep(NA, N_seq)
#
price_t_10_seq <- rep(NA, N_seq)
#
price_t_500_seq <- rep(NA, N_seq)
#
price_normal_seq <- rep(NA, N_seq)
#
Black_Scholes <- rep(NA, N_seq)

#
for(k in 1:N_seq){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T, 
                                               sigma = sigma,
                                               mu = mu,
                                               X_0 = X_0,
                                               K = K_seq[k], 
                                               r = r, 
                                               nu_hat = nu_hat, 
                                               delta = delta,
                                               N = N,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T,
                                          sigma = sigma,
                                          mu = mu,
                                          X_0 = X_0,
                                          K = K_seq[k], 
                                          r = r, 
                                          nu_hat = 3,
                                          delta = delta,
                                          N = N,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T, 
                                           sigma = sigma,
                                           mu = mu,
                                           X_0 = X_0,
                                           K = K_seq[k], 
                                           r = r, 
                                           nu_hat = 10, 
                                           delta = delta,
                                           N = N,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T, 
                                            sigma = sigma,
                                            mu = mu,
                                            X_0 = X_0,
                                            K = K_seq[k], 
                                            r = r, 
                                            nu_hat = 500,
                                            delta = delta,
                                            N = N,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T,
                                                  sigma = sigma,
                                                  mu = mu,
                                                  X_0 = X_0,
                                                  K = K_seq[k],
                                                  r = r,
                                                  delta = delta,
                                                  N = N)
  #
  Black_Scholes[k] <- BlackScholes(T = T,
                                   sigma = sigma,
                                   X_0 = X_0,
                                   K = K_seq[k],
                                   r = r)
}
#
#
par(mfrow = c(1,1))

#Legend
leg.txt.K <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(K_seq, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Strike Price"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Strike Price"),
        sub = TeX(sprintf("$X_0$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(X_0,4), signif(sigma,3), signif(r,3), signif(T,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("bottomleft", leg.txt.K, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)

#
# How option price depends on varying initial value
#
#######################

# Sequence of 20 possible initial values from 60 to 70
X_0_seq <- seq(from = 60,
               to = 70,
               length = N_seq)
#
price_t_nu_hat_seq <- rep(NA, N_seq)
#
price_t_3_seq <- rep(NA, N_seq)
#
price_t_10_seq <- rep(NA, N_seq)
#
price_t_500_seq <- rep(NA, N_seq)
#
price_normal_seq <- rep(NA, N_seq)
#
Black_Scholes <- rep(NA, N_seq)

#
for(k in 1:N_seq){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T, 
                                               sigma = sigma,
                                               mu = mu,
                                               X_0 = X_0_seq[k],
                                               K = K, 
                                               r = r, 
                                               nu_hat = nu_hat, 
                                               delta = delta,
                                               N = N,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T, 
                                          sigma = sigma,
                                          mu = mu,
                                          X_0 = X_0_seq[k],
                                          K = K, 
                                          r = r, 
                                          nu_hat = 3,
                                          delta = delta,
                                          N = N,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T, 
                                           sigma = sigma,
                                           mu = mu,
                                           X_0 = X_0_seq[k],
                                           K = K, 
                                           r = r, 
                                           nu_hat = 10, 
                                           delta = delta,
                                           N = N,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T, 
                                            sigma = sigma,
                                            mu = mu,
                                            X_0 = X_0_seq[k],
                                            K = K, 
                                            r = r, 
                                            nu_hat = 500, 
                                            delta = delta,
                                            N = N,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T,
                                                  sigma = sigma,
                                                  mu = mu,
                                                  X_0 = X_0_seq[k],
                                                  K = K,
                                                  r = r,
                                                  delta = delta,
                                                  N = N)
  #
  Black_Scholes[k] <- BlackScholes(T = T,
                                   sigma = sigma,
                                   X_0 = X_0_seq[k],
                                   K = K,
                                   r = r)
}
#
#
par(mfrow = c(1,1))

# Legend
leg.txt.X_0 <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(X_0_seq, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$X_0$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $X_0$"),
        sub = TeX(sprintf("$K$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(K,4), signif(sigma,3), signif(r,3), signif(T,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.X_0, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#dev.copy2pdf(file = "X_0.pdf",
#            width = 8,
#           height = 11)

#
# How option price depends on varying expiry time
#
#######################

# Sequence of 20 possible maturity values from a day to a year
T_seq <- seq(from = 1/252,
             to = 1,
             length = N_seq)
#
price_t_nu_hat_seq <- rep(NA, N_seq)
#
price_t_3_seq <- rep(NA, N_seq)
#
price_t_10_seq <- rep(NA, N_seq)
#
price_t_500_seq <- rep(NA, N_seq)
#
price_normal_seq <- rep(NA, N_seq)
#
Black_Scholes <- rep(NA, N_seq)

#
for(k in 1:N_seq){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_seq[k],
                                               sigma = sigma,
                                               mu = mu,
                                               X_0 = X_0,
                                               K = K, 
                                               r = r, 
                                               nu_hat = nu_hat, 
                                               delta = delta,
                                               N = N,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_seq[k],
                                          sigma = sigma,
                                          mu = mu,
                                          X_0 = X_0,
                                          K = K, 
                                          r = r, 
                                          nu_hat = 3, 
                                          delta = delta,
                                          N = N,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_seq[k], 
                                           sigma = sigma,
                                           mu = mu,
                                           X_0 = X_0,
                                           K = K, 
                                           r = r, 
                                           nu_hat = 10, 
                                           delta = delta,
                                           N = N,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_seq[k], 
                                            sigma = sigma,
                                            mu = mu,
                                            X_0 = X_0,
                                            K = K, 
                                            r = r, 
                                            nu_hat = 500, 
                                            delta = delta,
                                            N = N,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_seq[k],
                                                  sigma = sigma,
                                                  X_0 = X_0,
                                                  mu = mu,
                                                  K = K,
                                                  r = r,
                                                  delta = delta,
                                                  N = N)
  #
  Black_Scholes[k] <- BlackScholes(T = T_seq[k],
                                   sigma = sigma,
                                   X_0 = X_0,
                                   K = K,
                                   r = r)
}
#
#
par(mfrow = c(1,1))

# Legend
leg.txt.T <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(T_seq, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq,
              Black_Scholes),
        col = c("black", "green", "red", "turquoise", "blue","magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Time"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Time"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $r$ = %s", signif(X_0,4), signif(K,4), signif(sigma,3), signif(r,3))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.T, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)

#####################
# Appendix B

#
# JPM
#
#
# Interest Rate
#
#######################
# N_Seq <- 20
#
# T_2 <- 1
#
# K_2 <- JPM_Close_values[m_2] 
#
# sigma_2 <- sigma_hat_jpm
#
# X_0_2 <- JPM_Close_values[m_2-249]
#
# mu_2 <- mu_hat_jpm
#
# nu_hat_2 <- nu_hat_jpm # Degrees of freedom
#
r_seq_2 <- seq(from = 0.001,
               to = 0.05,
               length = 20)
#
N_2 <- 10000
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_t_normal_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = 1, # Years
                                               K = 138.63, 
                                               r = r_seq_2[k], 
                                               sigma = 0.2506884,
                                               X_0 = 97.11,
                                               mu = 0.1209486,
                                               nu_hat = 3.552891, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_2,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = 1, # Years
                                          K = 138.63, 
                                          r = r_seq_2[k], 
                                          sigma = 0.2506884,
                                          X_0 = 97.11,
                                          mu = 0.1209486,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_2,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = 1, # Years
                                           K = 138.63, 
                                           r = r_seq_2[k], 
                                           sigma = 0.2506884,
                                           X_0 = 97.11,
                                           mu = 0.1209486,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_2,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = 1, # Years
                                            K = 138.63, 
                                            r = r_seq_2[k], 
                                            sigma = 0.2506884,
                                            X_0 = 97.11,
                                            mu = 0.1209486,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_2,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = 1,
                                                  K = 138.63,
                                                  r = r_seq_2[k],
                                                  sigma = 0.2506884,
                                                  X_0 = 97.11,
                                                  mu = 0.1209486,
                                                  delta = 1 / 252,
                                                  N = N_2)
  #
  Black_Scholes[k] <- BlackScholes(T = 1,
                                   K = 138.63,
                                   r = r_seq_2[k],
                                   sigma = 0.2506884,
                                   X_0 = 97.11)
}
#
#
par(mfrow = c(1,1))

leg.txt.r <- c("t-distribution with 3 df", paste0("t estimated as ", signif(3.552891,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(r_seq_2, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Interest Rate"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Interest Rate"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $T$ = %s", signif(97.11,4), signif(138.63,4), signif(0.2506884,3), signif(1,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.r, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#
# Varying X_0
#
#######################
#
# r <- 0.005
#
X_0_seq_2 <- seq(from = 95,
                 to = 105,
                 length = 20)
#
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = 1, # Years
                                               K = 138.63, 
                                               r = 0.005, 
                                               sigma = 0.2506884,
                                               X_0 = X_0_seq_2[k],
                                               mu = 0.1209486,
                                               nu_hat = 3.552891, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_2,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = 1, # Years
                                          K = 138.63, 
                                          r = 0.005, 
                                          sigma = 0.2506884,
                                          X_0 = X_0_seq_2[k],
                                          mu = 0.1209486,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_2,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = 1, # Years
                                           K = 138.63, 
                                           r = 0.005, 
                                           sigma = 0.2506884,
                                           X_0 = X_0_seq_2[k],
                                           mu = 0.1209486,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_2,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = 1, # Years
                                            K = 138.63, 
                                            r = 0.005, 
                                            sigma = 0.2506884,
                                            X_0 = X_0_seq_2[k],
                                            mu = 0.1209486,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_2,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = 1,
                                                  K = 138.63,
                                                  r = 0.005,
                                                  sigma = 0.2506884,
                                                  X_0 = X_0_seq_2[k],
                                                  mu = 0.1209486,
                                                  delta = 1 / 252,
                                                  N = N_2)
  #
  Black_Scholes[k] <- BlackScholes(T = 1,
                                   K = 138.63,
                                   r = 0.005,
                                   sigma = 0.2506884,
                                   X_0 = X_0_seq_2[k])
}
#
#
par(mfrow = c(1,1))

leg.txt.X_0 <- c("t-distribution with 3 df", paste0("t estimated as ", signif(3.552891,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(X_0_seq_2, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$X_0$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $X_0$"),
        sub = TeX(sprintf("$K$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(138.63,4), signif(0.2506884,3), signif(0.005,3), signif(1,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.X_0, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#
# Varying sigma
#
#######################
#
#
sigma_seq_2 <- seq(from = 0.01,
                   to = 0.4,
                   length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = 1, # Years
                                               K = 138.63, 
                                               r = 0.005, 
                                               sigma = sigma_seq_2[k],
                                               X_0 = 97.11,
                                               mu = 0.1209486,
                                               nu_hat = 3.552891, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_2,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = 1, # Years
                                          K = 138.63, 
                                          r = 0.005, 
                                          sigma = sigma_seq_2[k],
                                          X_0 = 97.11,
                                          mu = 0.1209486,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_2,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = 1, # Years
                                           K = 138.63, 
                                           r = 0.005, 
                                           sigma = sigma_seq_2[k],
                                           X_0 = 97.11,
                                           mu = 0.1209486,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_2,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = 1, # Years
                                            K = 138.63, 
                                            r = 0.005, 
                                            sigma = sigma_seq_2[k],
                                            X_0 = 97.11,
                                            mu = 0.1209486,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_2,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = 1,
                                                  K = 138.63,
                                                  r = 0.005,
                                                  sigma = sigma_seq_2[k],
                                                  X_0 = 97.11,
                                                  mu = 0.1209486,
                                                  delta = 1 / 252,
                                                  N = N_2)
  #
  Black_Scholes[k] <- BlackScholes(T = 1,
                                   K = 138.63,
                                   r = 0.005,
                                   sigma = sigma_seq_2[k],
                                   X_0 = 97.11)
}
#
#
par(mfrow = c(1,1))

library(latex2exp)

leg.txt.sigma <- c("t-distribution with 3 df", paste0("t estimated as ", signif(3.552891,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes")
#
matplot(sigma_seq_2, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$\\sigma$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $\\sigma$"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $r$ = %s, $T$ = %s", signif(97.11,4), signif(138.63,4), signif(0.005,3), signif(1,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.sigma, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#
# Varying Strike Price
#
#######################
#
#
K_seq_2 <- seq(from = 97.11,
               to = 147.11,
               length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = 1, # Years
                                               K = K_seq_2[k], 
                                               r = 0.005, 
                                               sigma = 0.2506884,
                                               X_0 = 97.11,
                                               mu = 0.1209486,
                                               nu_hat = 3.552891, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_2,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = 1, # Years
                                          K = K_seq_2[k], 
                                          r = 0.005, 
                                          sigma = 0.2506884,
                                          X_0 = 97.11,
                                          mu = 0.1209486,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_2,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = 1, # Years
                                           K = K_seq_2[k], 
                                           r = 0.005, 
                                           sigma = 0.2506884,
                                           X_0 = 97.11,
                                           mu = 0.1209486,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_2,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = 1, # Years
                                            K = K_seq_2[k], 
                                            r = 0.005, 
                                            sigma = 0.2506884,
                                            X_0 = 97.11,
                                            mu = 0.1209486,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_2,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = 1,
                                                  K = K_seq_2[k],
                                                  r = 0.005,
                                                  sigma = 0.2506884,
                                                  X_0 = 97.11,
                                                  mu = 0.1209486,
                                                  delta = 1 / 252,
                                                  N = N_2)
  #
  Black_Scholes[k] <- BlackScholes(T = 1,
                                   K = K_seq_2[k],
                                   r = 0.005,
                                   sigma = 0.2506884,
                                   X_0 = 97.11)
}
#
#
par(mfrow = c(1,1))

leg.txt.K <- c("t-distribution with 3 df", paste0("t estimated as ", signif(3.552891,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(K_seq_2, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Strike Price"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Strike Price"),
        sub = TeX(sprintf("$X_0$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(97.11,4), signif(0.2506884,3), signif(0.005,3), signif(1,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("bottomleft", leg.txt.K, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)

#
# Varying Expiry Time
#
#######################
#
#
T_seq_2 <- seq(from = 0.01,
               to = 1,
               length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_seq_2[k], # Years
                                               K = 138.63, 
                                               r = 0.005, 
                                               sigma = 0.2506884,
                                               X_0 = 97.11,
                                               mu = 0.1209486,
                                               nu_hat = 3.552891, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_2,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_seq_2[k], # Years
                                          K = 138.63, 
                                          r = 0.005, 
                                          sigma = 0.2506884,
                                          X_0 = 97.11,
                                          mu = 0.1209486,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_2,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_seq_2[k], # Years
                                           K = 138.63, 
                                           r = 0.005, 
                                           sigma = 0.2506884,
                                           X_0 = 97.11,
                                           mu = 0.1209486,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_2,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_seq_2[k], # Years
                                            K = 138.63, 
                                            r = 0.005, 
                                            sigma = 0.2506884,
                                            X_0 = 97.11,
                                            mu = 0.1209486,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_2,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_seq_2[k],
                                                  K = 138.63,
                                                  r = 0.005,
                                                  sigma = 0.2506884,
                                                  X_0 = 97.11,
                                                  mu = 0.1209486,
                                                  delta = 1/252,
                                                  N = N_2)
  #
  Black_Scholes[k] <- BlackScholes(T = T_seq_2[k],
                                   K = 138.63,
                                   r = 0.005,
                                   sigma = 0.2506884,
                                   X_0 = 97.11)
}
#
#
par(mfrow = c(1,1))

leg.txt.T <- c("t-distribution with 3 df", paste0("t estimated as ", signif(3.552891,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(T_seq_2, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq,
              Black_Scholes),
        col = c("black", "green", "red", "turquoise", "blue","magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Time"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Time"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $r$ = %s", signif(97.11,4), signif(138.63,4), signif(0.2506884,3), signif(0.005,3))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.T, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)


################################
#
# Atlas Copco
#
#
# Interest Rate
#
#######################
#
#
T_3 <- 1
#
K_3 <- A_C_Close_values[m_3+1] 
#
sigma_3 <- sigma_hat_A_C
#
X_0_3 <- A_C_Close_values[m_3-245]
#
mu_3 <- mu_hat_A_C
#
nu_hat_3 <- nu_hat_A_C # Degrees of freedom
#
r_seq_3 <- seq(from = 0.001,
               to = 0.05,
               length = 20)
#
N_3 <- 10000
#
r_constant <- 0.005
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_t_normal_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_3, # Years
                                               K = K_3, 
                                               r = r_seq_3[k], 
                                               sigma = sigma_3,
                                               X_0 = X_0_3,
                                               mu = mu_3,
                                               nu_hat = nu_hat_3, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_3,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_3, # Years
                                          K = K_3, 
                                          r = r_seq_3[k], 
                                          sigma = sigma_3,
                                          X_0 = X_0_3,
                                          mu = mu_3,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_3,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_3, # Years
                                           K = K_3, 
                                           r = r_seq_3[k], 
                                           sigma = sigma_3,
                                           X_0 = X_0_3,
                                           mu = mu_3,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_3,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_3, # Years
                                            K = K_3, 
                                            r = r_seq_3[k], 
                                            sigma = sigma_3,
                                            X_0 = X_0_3,
                                            mu = mu_3,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_3,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_3,
                                                  K = K_3,
                                                  r = r_seq_3[k],
                                                  sigma = sigma_3,
                                                  X_0 = X_0_3,
                                                  mu = mu_3,
                                                  delta = 1 / 252,
                                                  N = N_3)
  #
  Black_Scholes[k] <- BlackScholes(T = T_3,
                                   K = K_3,
                                   r = r_seq_3[k],
                                   sigma = sigma_3,
                                   X_0 = X_0_3)
}
#
#
par(mfrow = c(1,1))

leg.txt.r <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat_3,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(r_seq_3, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Interest Rate"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Interest Rate"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $T$ = %s", signif(X_0_3,4), signif(K_3,4), signif(sigma_3,3), signif(T_3,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.r, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

# Varying X_0
#
#######################
#
#
X_0_seq_3 <- seq(from = 45,
                 to = 65,
                 length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_3, # Years
                                               K = K_3, 
                                               r = 0.005, 
                                               sigma = sigma_3,
                                               X_0 = X_0_seq_3[k],
                                               mu = mu_3,
                                               nu_hat = nu_hat_3, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_3,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_3, # Years
                                          K = K_3, 
                                          r = 0.005, 
                                          sigma = sigma_3,
                                          X_0 = X_0_seq_3[k],
                                          mu = mu_3,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_3,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_3, # Years
                                           K = K_3, 
                                           r = 0.005, 
                                           sigma = sigma_3,
                                           X_0 = X_0_seq_3[k],
                                           mu = mu_3,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_3,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_3, # Years
                                            K = K_3, 
                                            r = 0.005, 
                                            sigma = sigma_3,
                                            X_0 = X_0_seq_3[k],
                                            mu = mu_3,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_3,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_3,
                                                  K = K_3,
                                                  r = 0.005,
                                                  sigma = sigma_3,
                                                  X_0 = X_0_seq_3[k],
                                                  mu = mu_3,
                                                  delta = 1 / 252,
                                                  N = N_3)
  #
  Black_Scholes[k] <- BlackScholes(T = T_3,
                                   K = K_3,
                                   r = 0.005,
                                   sigma = sigma_3,
                                   X_0 = X_0_seq_3[k])
}
#
#
par(mfrow = c(1,1))

leg.txt.X_0 <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat_3,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(X_0_seq_3, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$X_0$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $X_0$"),
        sub = TeX(sprintf("$K$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(K_3,4), signif(sigma_3,3), signif(0.005,3), signif(T_3,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.X_0, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#
# Varying sigma
#
#######################
#
#
sigma_seq_3 <- seq(from = 0.01,
                   to = 0.4,
                   length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_3, # Years
                                               K = K_3, 
                                               r = 0.005, 
                                               sigma = sigma_seq_3[k],
                                               X_0 = X_0_3,
                                               mu = mu_3,
                                               nu_hat = nu_hat_3, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_3,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_3, # Years
                                          K = K_3, 
                                          r = 0.005, 
                                          sigma = sigma_seq_3[k],
                                          X_0 = X_0_3,
                                          mu = mu_3,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_3,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_3, # Years
                                           K = K_3, 
                                           r = 0.005, 
                                           sigma = sigma_seq_3[k],
                                           X_0 = X_0_3,
                                           mu = mu_3,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_3,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_3, # Years
                                            K = K_3, 
                                            r = 0.005, 
                                            sigma = sigma_seq_3[k],
                                            X_0 = X_0_3,
                                            mu = mu_3,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_3,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_3,
                                                  K = K_3,
                                                  r = 0.005,
                                                  sigma = sigma_seq_3[k],
                                                  X_0 = X_0_3,
                                                  mu = mu_3,
                                                  delta = 1 / 252,
                                                  N = N_3)
  #
  Black_Scholes[k] <- BlackScholes(T = T_3,
                                   K = K_3,
                                   r = 0.005,
                                   sigma = sigma_seq_3[k],
                                   X_0 = X_0_3)
}
#
#
par(mfrow = c(1,1))

library(latex2exp)

leg.txt.sigma <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat_3,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes")
#
matplot(sigma_seq_3, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("$\\sigma$"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on $\\sigma$"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $r$ = %s, $T$ = %s", signif(X_0_3,4), signif(K_3,4), signif(0.005,3), signif(T_3,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.sigma, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue", "magenta"),
       cex = 1)

#
# Varying Strike Price
#
#######################
#
#
K_seq_3 <- seq(from = X_0_3,
               to = X_0_3+50,
               length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_3, # Years
                                               K = K_seq_3[k], 
                                               r = 0.005, 
                                               sigma = sigma_3,
                                               X_0 = X_0_3,
                                               mu = mu_3,
                                               nu_hat = nu_hat_3, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_3,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_3, # Years
                                          K = K_seq_3[k], 
                                          r = 0.005, 
                                          sigma = sigma_3,
                                          X_0 = X_0_3,
                                          mu = mu_3,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_3,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_3, # Years
                                           K = K_seq_3[k], 
                                           r = 0.005, 
                                           sigma = sigma_3,
                                           X_0 = X_0_3,
                                           mu = mu_3,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_3,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_3, # Years
                                            K = K_seq_3[k], 
                                            r = 0.005, 
                                            sigma = sigma_3,
                                            X_0 = X_0_3,
                                            mu = mu_3,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_3,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_3,
                                                  K = K_seq_3[k],
                                                  r = 0.005,
                                                  sigma = sigma_3,
                                                  X_0 = X_0_3,
                                                  mu = mu_3,
                                                  delta = 1 / 252,
                                                  N = N_3)
  #
  Black_Scholes[k] <- BlackScholes(T = T_3,
                                   K = K_seq_3[k],
                                   r = 0.005,
                                   sigma = sigma_3,
                                   X_0 = X_0_3)
}
#
#
par(mfrow = c(1,1))

leg.txt.K <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat_3,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(K_seq_3, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq),
        col = c("black", "green", "red", "turquoise", "blue", "magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Strike Price"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Strike Price"),
        sub = TeX(sprintf("$X_0$ = %s, $sigma$ = %s, $r$ = %s, $T$ = %s", signif(X_0_3,4), signif(sigma_3,3), signif(0.005,3), signif(T_3,1))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("bottomleft", leg.txt.K, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)

#
# Varying Expiry Time
#
#######################
#
#
T_seq_3 <- seq(from = 0.01,
               to = 1,
               length = 20)
#
price_t_nu_hat_seq <- rep(NA, 20)
#
price_t_3_seq <- rep(NA, 20)
#
price_t_10_seq <- rep(NA, 20)
#
price_t_500_seq <- rep(NA, 20)
#
price_normal_seq <- rep(NA, 20)
#
Black_Scholes <- rep(NA, 20)
#
#
for(k in 1:20){
  #
  print(k)
  #
  price_t_nu_hat_seq[k] <- option_t_increments(T = T_seq_3[k], # Years
                                               K = K_3, 
                                               r = 0.005, 
                                               sigma = sigma_3,
                                               X_0 = X_0_3,
                                               mu = mu_3,
                                               nu_hat = nu_hat_3, #Degrees of freedom
                                               delta = 1 / 252,
                                               N = N_3,
                                               do_trim = TRUE, 
                                               trim = 0.05)
  #
  price_t_3_seq[k] <- option_t_increments(T = T_seq_3[k], # Years
                                          K = K_3, 
                                          r = 0.005, 
                                          sigma = sigma_3,
                                          X_0 = X_0_3,
                                          mu = mu_3,
                                          nu_hat = 3, # Degrees of freedom
                                          delta = 1 / 252,
                                          N = N_3,
                                          do_trim = TRUE, 
                                          trim = 0.05)
  #
  price_t_10_seq[k] <- option_t_increments(T = T_seq_3[k], # Years
                                           K = K_3, 
                                           r = 0.005, 
                                           sigma = sigma_3,
                                           X_0 = X_0_3,
                                           mu = mu_3,
                                           nu_hat = 10, # Degrees of freedom
                                           delta = 1 / 252,
                                           N = N_3,
                                           do_trim = TRUE, 
                                           trim = 0.05)
  #
  price_t_500_seq[k] <- option_t_increments(T = T_seq_3[k], # Years
                                            K = K_3, 
                                            r = 0.005, 
                                            sigma = sigma_3,
                                            X_0 = X_0_3,
                                            mu = mu_3,
                                            nu_hat = 500, # Degrees of freedom
                                            delta = 1 / 252,
                                            N = N_3,
                                            do_trim = TRUE, 
                                            trim = 0.05)
  #
  price_normal_seq[k] <- option_normal_increments(T = T_seq_3[k],
                                                  K = K_3,
                                                  r = 0.005,
                                                  sigma = sigma_3,
                                                  X_0 = X_0_3,
                                                  mu = mu_3,
                                                  delta = 1 / 252,
                                                  N = N_3)
  #
  Black_Scholes[k] <- BlackScholes(T = T_seq_3[k],
                                   K = K_3,
                                   r = 0.005,
                                   sigma = sigma_3,
                                   X_0 = X_0_3)
}
#
#
par(mfrow = c(1,1))

leg.txt.T <- c("t-distribution with 3 df", paste0("t estimated as ", signif(nu_hat_3,2), " df"),  "t-distribution with 10 df", "t-distribution with 500 df", "Normal", "Black-Scholes Model")
#
matplot(T_seq_3, 
        cbind(price_t_3_seq,
              price_t_nu_hat_seq,
              price_t_10_seq,
              price_t_500_seq,
              price_normal_seq,
              Black_Scholes),
        col = c("black", "green", "red", "turquoise", "blue","magenta"),
        type = "b",
        lty = 1,
        lwd = 1.5,
        pch = 4,
        xlab = TeX("Time"),
        ylab = "Option Price",
        main = TeX("The Dependence of Option Price on Time"),
        sub = TeX(sprintf("$X_0$ = %s, $K$ = %s, $sigma$ = %s, $r$ = %s", signif(X_0_3,4), signif(K_3,4), signif(sigma_3,3), signif(0.005,3))),
        cex.lab = 1.5, # Expand the labels
        cex.axis = 1, # Expand the numbers
        cex.main = 1.5, # Expand the main title
        cex.sub = 0.8) # Expand the subtitle
legend("topleft", leg.txt.T, pch=4, title= "Distribution", 
       lty = 1,  
       lwd = 1.5,
       col = c("black", "green", "red", "turquoise","blue","magenta"),
       cex = 1)
