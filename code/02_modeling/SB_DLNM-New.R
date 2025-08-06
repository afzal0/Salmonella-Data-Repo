###############################################################################
# SINGLE R SCRIPT: MODELS 3 & 4 WITH NIMBLE + iCAR
# Using for-loops for crossbasis sums (no partial indexing in inprod())
###############################################################################

# 1) LOAD LIBRARIES & DATA
library(nimble)
library(dplyr)
library(splines)
library(dlnm)

# Load your data
load("model_input/data_casecrossover.RData") # => data_cco
load("model_input/data_timeseries.RData")    # => data_ts
load("model_input/spatial_adjacency.RData")  # => adj_vec, wts_vec, num_vec

# Ensure no row names
rownames(data_cco) <- NULL
rownames(data_ts)  <- NULL

# Fix any duplicate column names, since dplyr >= 1.1.0 forbids mutate on duplicates
if(any(duplicated(names(data_cco)))) {
  names(data_cco) <- make.unique(names(data_cco))
}
if(any(duplicated(names(data_ts)))) {
  names(data_ts) <- make.unique(names(data_ts))
}

# 2) CENTER EXPOSURES
data_cco <- data_cco %>%
  mutate(
    Tmean_centered     = Tmean - mean(Tmean, na.rm=TRUE),
    Rainfall_centered  = Rainfall - mean(Rainfall, na.rm=TRUE),
    Floodwave_centered = Floodwave - mean(Floodwave, na.rm=TRUE),
    Humidity_centered  = Humidity - mean(Humidity, na.rm=TRUE)
  )

data_ts <- data_ts %>%
  mutate(
    Tmean_centered     = Tmean - mean(Tmean, na.rm=TRUE),
    Rainfall_centered  = Rainfall - mean(Rainfall, na.rm=TRUE),
    Floodwave_centered = Floodwave - mean(Floodwave, na.rm=TRUE),
    Humidity_centered  = Humidity - mean(Humidity, na.rm=TRUE)
  )

# 3) CREATE CROSSBASIS
max_lag   <- 2
var_knots <- c(0.10, 0.50, 0.90)

## For Model 3 (Case-Crossover)
cb_cco_Tmean <- crossbasis(data_cco$Tmean_centered, lag=max_lag,
                           argvar=list(fun="ns", knots=quantile(data_cco$Tmean_centered, var_knots, na.rm=TRUE)),
                           arglag=list(fun="ns"))
cb_cco_Rain <- crossbasis(data_cco$Rainfall_centered, lag=max_lag,
                          argvar=list(fun="ns", knots=quantile(data_cco$Rainfall_centered, var_knots, na.rm=TRUE)),
                          arglag=list(fun="ns"))
cb_cco_Flood <- crossbasis(data_cco$Floodwave_centered, lag=max_lag,
                           argvar=list(fun="ns", knots=quantile(data_cco$Floodwave_centered, var_knots, na.rm=TRUE)),
                           arglag=list(fun="ns"))
cb_cco_Humid <- crossbasis(data_cco$Humidity_centered, lag=max_lag,
                           argvar=list(fun="ns", knots=quantile(data_cco$Humidity_centered, var_knots, na.rm=TRUE)),
                           arglag=list(fun="ns"))

my_wtemp_cco     <- as.matrix(cb_cco_Tmean)
my_wrain_cco     <- as.matrix(cb_cco_Rain)
my_wflood_cco    <- as.matrix(cb_cco_Flood)
my_whumidity_cco <- as.matrix(cb_cco_Humid)

## For Model 4 (Time-Series)
cb_ts_Tmean <- crossbasis(data_ts$Tmean_centered, lag=max_lag,
                          argvar=list(fun="ns", knots=quantile(data_ts$Tmean_centered, var_knots, na.rm=TRUE)),
                          arglag=list(fun="ns"))
cb_ts_Rain <- crossbasis(data_ts$Rainfall_centered, lag=max_lag,
                         argvar=list(fun="ns", knots=quantile(data_ts$Rainfall_centered, var_knots, na.rm=TRUE)),
                         arglag=list(fun="ns"))
cb_ts_Flood <- crossbasis(data_ts$Floodwave_centered, lag=max_lag,
                          argvar=list(fun="ns", knots=quantile(data_ts$Floodwave_centered, var_knots, na.rm=TRUE)),
                          arglag=list(fun="ns"))
cb_ts_Humid <- crossbasis(data_ts$Humidity_centered, lag=max_lag,
                          argvar=list(fun="ns", knots=quantile(data_ts$Humidity_centered, var_knots, na.rm=TRUE)),
                          arglag=list(fun="ns"))

my_wtemp_ts     <- as.matrix(cb_ts_Tmean)
my_wrain_ts     <- as.matrix(cb_ts_Rain)
my_wflood_ts    <- as.matrix(cb_ts_Flood)
my_whumidity_ts <- as.matrix(cb_ts_Humid)

# 4) CREATE OUTPUT DIRECTORY IF NEEDED
dir.create("model_output", showWarnings=FALSE)

# 5) MCMC PARAMETERS
n.chains <- 3
n.iter   <- 5000

###############################################################################
# MODEL 3: CASE-CROSSOVER + iCAR (NIMBLE)
###############################################################################
# 1) Prepare data & constants
model3_data <- list(
  y = as.numeric(data_cco$Salmonella)
)

code_model3 <- nimbleCode({
  
  #############################################################################
  # 1) Hyperpriors for crossbasis - unchanged
  #############################################################################
  
  # Tmean
  for(j in 1:n_coef_temp) {
    mu_temp[j]    ~ dnorm(0, 0.001)
    sigma_temp[j] ~ dunif(0, 10)
  }
  tau_temp ~ dgamma(0.1, 0.1)
  
  # Rain
  for(j in 1:n_coef_rain) {
    mu_rain[j]    ~ dnorm(0, 0.001)
    sigma_rain[j] ~ dunif(0, 10)
  }
  tau_rain ~ dgamma(0.1, 0.1)
  
  # Flood
  for(j in 1:n_coef_flood) {
    mu_flood[j]    ~ dnorm(0, 0.001)
    sigma_flood[j] ~ dunif(0, 10)
  }
  tau_flood ~ dgamma(0.1, 0.1)
  
  # Humidity
  for(j in 1:n_coef_humidity) {
    mu_humidity[j]    ~ dnorm(0, 0.001)
    sigma_humidity[j] ~ dunif(0, 10)
  }
  tau_humidity ~ dgamma(0.1, 0.1)
  
  #############################################################################
  # 2) Stratum intercepts - unchanged
  #############################################################################
  for(s in 1:n_strata) {
    alpha[s] ~ dnorm(0, 0.001)
  }
  
  #############################################################################
  # 3) iCAR Spatial - unchanged
  #############################################################################
  tau_phi ~ dgamma(0.1, 0.1)
  phi[1:n_reg] ~ dcar_normal(adj[1:n_adj], weights[1:n_adj], num[1:n_reg],
                             tau_phi, zero_mean=1)
  
  #############################################################################
  # 4) Region-level random effects for crossbasis - unchanged
  #############################################################################
  for(r in 1:n_reg) {
    # Tmean
    for(j in 1:n_coef_temp) {
      theta_temp[r,j] ~ dnorm(0, tau_temp)
      beta_temp[r,j]  <- mu_temp[j] + sigma_temp[j] * theta_temp[r,j]
    }
    # Rain
    for(j in 1:n_coef_rain) {
      theta_rain[r,j] ~ dnorm(0, tau_rain)
      beta_rain[r,j]  <- mu_rain[j] + sigma_rain[j] * theta_rain[r,j]
    }
    # Flood
    for(j in 1:n_coef_flood) {
      theta_flood[r,j] ~ dnorm(0, tau_flood)
      beta_flood[r,j]  <- mu_flood[j] + sigma_flood[j] * theta_flood[r,j]
    }
    # Humidity
    for(j in 1:n_coef_humidity) {
      theta_humidity[r,j] ~ dnorm(0, tau_humidity)
      beta_humidity[r,j]  <- mu_humidity[j] + sigma_humidity[j] * theta_humidity[r,j]
    }
  }
  
  #############################################################################
  # 5) Likelihood using inprod() instead of for-loops
  #############################################################################
  for(i in 1:n_data) {
    y[i] ~ dpois(lambda[i])
    
    # Calculate crossbasis contributions using inprod()
    for(j in 1:n_coef_temp) {
      beta_temp_i[i,j] <- beta_temp[region_index[i], j]
    }
    temp_contrib[i] <- inprod(beta_temp_i[i,1:n_coef_temp], w_temp[i,1:n_coef_temp])
    
    for(j in 1:n_coef_rain) {
      beta_rain_i[i,j] <- beta_rain[region_index[i], j]
    }
    rain_contrib[i] <- inprod(beta_rain_i[i,1:n_coef_rain], w_rain[i,1:n_coef_rain])
    
    for(j in 1:n_coef_flood) {
      beta_flood_i[i,j] <- beta_flood[region_index[i], j]
    }
    flood_contrib[i] <- inprod(beta_flood_i[i,1:n_coef_flood], w_flood[i,1:n_coef_flood])
    
    for(j in 1:n_coef_humidity) {
      beta_humidity_i[i,j] <- beta_humidity[region_index[i], j]
    }
    humid_contrib[i] <- inprod(beta_humidity_i[i,1:n_coef_humidity], w_humidity[i,1:n_coef_humidity])
    
    log(lambda[i]) <- alpha[strata[i]] +
      temp_contrib[i] + rain_contrib[i] +
      flood_contrib[i] + humid_contrib[i] +
      phi[region_index[i]]
    
    dev_i[i] <- -2 * log(dpois(y[i], lambda[i]) + 1e-15)
  }
  
  # Summation of deviance
  D <- sum(dev_i[1:n_data])
})

# Update constants to include these new indexed variables
model3_constants <- list(
  n_data   = nrow(data_cco),
  n_strata = length(unique(data_cco$strata)),
  n_reg    = length(unique(data_cco$LHD_code)),
  
  # Indices
  strata       = as.numeric(factor(data_cco$strata)),
  region_index = as.numeric(data_cco$LHD_code),
  
  # Crossbasis
  w_temp     = my_wtemp_cco,
  w_rain     = my_wrain_cco,
  w_flood    = my_wflood_cco,
  w_humidity = my_whumidity_cco,
  
  n_coef_temp     = ncol(my_wtemp_cco),
  n_coef_rain     = ncol(my_wrain_cco),
  n_coef_flood    = ncol(my_wflood_cco),
  n_coef_humidity = ncol(my_whumidity_cco),
  
  # iCAR adjacency
  adj     = adj_vec,
  weights = wts_vec,
  num     = num_vec,
  n_adj   = length(adj_vec)
)

# 3) Initial values
inits_model3 <- function() {
  list(
    alpha = rep(0, model3_constants$n_strata),
    
    mu_temp     = rep(0, model3_constants$n_coef_temp),
    sigma_temp  = rep(1, model3_constants$n_coef_temp),
    tau_temp    = 1,
    
    mu_rain     = rep(0, model3_constants$n_coef_rain),
    sigma_rain  = rep(1, model3_constants$n_coef_rain),
    tau_rain    = 1,
    
    mu_flood    = rep(0, model3_constants$n_coef_flood),
    sigma_flood = rep(1, model3_constants$n_coef_flood),
    tau_flood   = 1,
    
    mu_humidity    = rep(0, model3_constants$n_coef_humidity),
    sigma_humidity = rep(1, model3_constants$n_coef_humidity),
    tau_humidity   = 1,
    
    tau_phi = 1,
    phi     = rep(0, model3_constants$n_reg),
    
    theta_temp     = matrix(0, model3_constants$n_reg, model3_constants$n_coef_temp),
    theta_rain     = matrix(0, model3_constants$n_reg, model3_constants$n_coef_rain),
    theta_flood    = matrix(0, model3_constants$n_reg, model3_constants$n_coef_flood),
    theta_humidity = matrix(0, model3_constants$n_reg, model3_constants$n_coef_humidity)
  )
}

# 4) Build & compile model
Rmodel3 <- nimbleModel(code=code_model3,
                       data=model3_data,
                       constants=model3_constants,
                       inits=inits_model3())
Cmodel3 <- compileNimble(Rmodel3)

# 5) Configure MCMC
spec3 <- buildMCMC(Rmodel3,
                   monitors=c("alpha","beta_temp","beta_rain","beta_flood","beta_humidity","phi","D"))
Cspec3 <- compileNimble(spec3, project=Rmodel3)

# 6) Run MCMC
cat("\n--- Running Model 3 (Case-Crossover + iCAR) ---\n")
Cspec3$run(n.iter)
model3_samples <- as.matrix(Cspec3$mvSamples)
save(model3_samples, file="model_output/model3_results_nimble.RData")
cat("Model 3 completed and saved.\n")

###############################################################################
# MODEL 4: TIME-SERIES + iCAR (NIMBLE)
###############################################################################
# 1) Trend & Seasonal


trend_basis <- ns(1:nrow(data_ts), df=4)

# Example: a natural spline basis for month (for seasonal pattern)
month_basis <- ns(1:12, df=3)

# Expand that monthly basis to each row in the time series
seas_basis  <- matrix(NA, nrow=nrow(data_ts), ncol=ncol(month_basis))
for(i in seq_len(nrow(data_ts))) {
  seas_basis[i, ] <- month_basis[data_ts$month[i], ]
}

# (3) Prepare data and constants for NIMBLE
model4_data <- list(
  y = as.numeric(data_ts$Salmonella)
)

# The model code, with a single-line definition for `season_val[i]`
# to avoid multiple definitions:
code_model4 <- nimbleCode({
  
  ############################
  # Crossbasis priors
  ############################
  # Tmean
  for(j in 1:n_coef_temp) {
    mu_temp[j]    ~ dnorm(0, 0.001)
    sigma_temp[j] ~ dunif(0, 10)
  }
  tau_temp ~ dgamma(0.1, 0.1)
  
  # Rain
  for(j in 1:n_coef_rain) {
    mu_rain[j]    ~ dnorm(0, 0.001)
    sigma_rain[j] ~ dunif(0, 10)
  }
  tau_rain ~ dgamma(0.1, 0.1)
  
  # Flood
  for(j in 1:n_coef_flood) {
    mu_flood[j]    ~ dnorm(0, 0.001)
    sigma_flood[j] ~ dunif(0, 10)
  }
  tau_flood ~ dgamma(0.1, 0.1)
  
  # Humidity
  for(j in 1:n_coef_humidity) {
    mu_humidity[j]    ~ dnorm(0, 0.001)
    sigma_humidity[j] ~ dunif(0, 10)
  }
  tau_humidity ~ dgamma(0.1, 0.1)
  
  ############################
  # Region-level intercepts by month
  ############################
  for(r in 1:n_reg) {
    for(m in 1:12) {
      alpha[r,m] ~ dnorm(0, 0.001)
    }
  }
  
  ############################
  # iCAR
  ############################
  tau_phi ~ dgamma(0.1, 0.1)
  phi[1:n_reg] ~ dcar_normal(
    adj[1:n_adj],
    weights[1:n_adj],
    num[1:n_reg],
    tau_phi,
    zero_mean = 1
  )
  
  ############################
  # Crossbasis random effects
  ############################
  for(r in 1:n_reg) {
    # Tmean
    for(j in 1:n_coef_temp) {
      theta_temp[r,j] ~ dnorm(0, tau_temp)
      beta_temp[r,j]  <- mu_temp[j] + sigma_temp[j] * theta_temp[r,j]
    }
    # Rain
    for(j in 1:n_coef_rain) {
      theta_rain[r,j] ~ dnorm(0, tau_rain)
      beta_rain[r,j]  <- mu_rain[j] + sigma_rain[j] * theta_rain[r,j]
    }
    # Flood
    for(j in 1:n_coef_flood) {
      theta_flood[r,j] ~ dnorm(0, tau_flood)
      beta_flood[r,j]  <- mu_flood[j] + sigma_flood[j] * theta_flood[r,j]
    }
    # Humidity
    for(j in 1:n_coef_humidity) {
      theta_humidity[r,j] ~ dnorm(0, tau_humidity)
      beta_humidity[r,j]  <- mu_humidity[j] + sigma_humidity[j] * theta_humidity[r,j]
    }
  }
  
  ############################
  # Trend & seasonal random effects
  ############################
  for(r in 1:n_reg) {
    # Trend
    for(tt in 1:n_trend) {
      gamma[r,tt] ~ dnorm(0, 0.001)
    }
    # Seasonal differences by year
    for(yy in 1:n_year) {
      for(ssi in 1:n_seas) {
        delta[r,yy,ssi] ~ dnorm(0, 0.001)
      }
    }
  }
  
  ############################
  # Likelihood using inprod()
  ############################
  for(i in 1:n_data) {
    y[i] ~ dpois(lambda[i])
    
    # Crossbasis: Tmean
    for(j in 1:n_coef_temp) {
      beta_temp_i[i,j] <- beta_temp[region_index[i], j]
    }
    temp_contrib[i] <- inprod(beta_temp_i[i,1:n_coef_temp],
                              w_temp[i,1:n_coef_temp])
    
    # Crossbasis: Rain
    for(j in 1:n_coef_rain) {
      beta_rain_i[i,j] <- beta_rain[region_index[i], j]
    }
    rain_contrib[i] <- inprod(beta_rain_i[i,1:n_coef_rain],
                              w_rain[i,1:n_coef_rain])
    
    # Crossbasis: Flood
    for(j in 1:n_coef_flood) {
      beta_flood_i[i,j] <- beta_flood[region_index[i], j]
    }
    flood_contrib[i] <- inprod(beta_flood_i[i,1:n_coef_flood],
                               w_flood[i,1:n_coef_flood])
    
    # Crossbasis: Humidity
    for(j in 1:n_coef_humidity) {
      beta_humidity_i[i,j] <- beta_humidity[region_index[i], j]
    }
    humid_contrib[i] <- inprod(beta_humidity_i[i,1:n_coef_humidity],
                               w_humidity[i,1:n_coef_humidity])
    
    # Trend contribution
    for(j in 1:n_trend) {
      gamma_i[i,j] <- gamma[region_index[i], j]
    }
    trend_val[i] <- inprod(gamma_i[i,1:n_trend],
                           t[i,1:n_trend])
    
    # Seasonal contribution in one statement
    season_val[i] <- inprod(delta[region_index[i], year[i], 1:n_seas],
                            s[i, 1:n_seas])
    
    # Linear predictor
    log(lambda[i]) <-
      alpha[region_index[i], month[i]] +
      temp_contrib[i] + rain_contrib[i] +
      flood_contrib[i] + humid_contrib[i] +
      trend_val[i] + season_val[i] +
      phi[region_index[i]]
    
    # Deviance
    dev_i[i] <- -2 * log(dpois(y[i], lambda[i]) + 1e-15)
  }
  # Summation of deviance
  D <- sum(dev_i[1:n_data])
})

# (4) Model constants
model4_constants <- list(
  n_data = nrow(data_ts),
  n_reg  = length(unique(data_ts$LHD_code)),
  
  region_index = as.numeric(data_ts$LHD_code),
  month        = as.numeric(data_ts$month),
  year         = as.numeric(data_ts$year - min(data_ts$year) + 1),
  
  w_temp     = my_wtemp_ts,
  w_rain     = my_wrain_ts,
  w_flood    = my_wflood_ts,
  w_humidity = my_whumidity_ts,
  
  n_coef_temp     = ncol(my_wtemp_ts),
  n_coef_rain     = ncol(my_wrain_ts),
  n_coef_flood    = ncol(my_wflood_ts),
  n_coef_humidity = ncol(my_whumidity_ts),
  
  t       = trend_basis,
  s       = seas_basis,
  n_trend = ncol(trend_basis),
  n_seas  = ncol(seas_basis),
  n_year  = length(unique(data_ts$year)),
  
  # iCAR adjacency
  adj     = adj_vec,
  weights = wts_vec,
  num     = num_vec,
  n_adj   = length(adj_vec)
)

# (5) Initial values function
inits_model4 <- function() {
  list(
    mu_temp     = rep(0, model4_constants$n_coef_temp),
    sigma_temp  = rep(1, model4_constants$n_coef_temp),
    tau_temp    = 1,
    
    mu_rain     = rep(0, model4_constants$n_coef_rain),
    sigma_rain  = rep(1, model4_constants$n_coef_rain),
    tau_rain    = 1,
    
    mu_flood    = rep(0, model4_constants$n_coef_flood),
    sigma_flood = rep(1, model4_constants$n_coef_flood),
    tau_flood   = 1,
    
    mu_humidity    = rep(0, model4_constants$n_coef_humidity),
    sigma_humidity = rep(1, model4_constants$n_coef_humidity),
    tau_humidity   = 1,
    
    alpha = matrix(0, model4_constants$n_reg, 12),
    
    tau_phi = 1,
    phi     = rep(0, model4_constants$n_reg),
    
    theta_temp     = matrix(0, model4_constants$n_reg, model4_constants$n_coef_temp),
    theta_rain     = matrix(0, model4_constants$n_reg, model4_constants$n_coef_rain),
    theta_flood    = matrix(0, model4_constants$n_reg, model4_constants$n_coef_flood),
    theta_humidity = matrix(0, model4_constants$n_reg, model4_constants$n_coef_humidity),
    
    gamma = matrix(0, model4_constants$n_reg, model4_constants$n_trend),
    delta = array(0, dim = c(model4_constants$n_reg,
                             model4_constants$n_year,
                             model4_constants$n_seas))
  )
}

# (6) Build & compile the model
Rmodel4 <- nimbleModel(
  code      = code_model4,
  data      = model4_data,
  constants = model4_constants,
  inits     = inits_model4()
)
Cmodel4 <- compileNimble(Rmodel4)

# (7) Configure and compile the MCMC
spec4 <- buildMCMC(
  Rmodel4,
  monitors = c("alpha","beta_temp","beta_rain","beta_flood","beta_humidity",
               "gamma","delta","phi","D")
)
Cspec4 <- compileNimble(spec4, project = Rmodel4)

# (8) Run MCMC
n.chains <- 3
n.iter   <- 5000  # Example
cat("\n--- Running Model 4 (Time-Series + iCAR) ---\n")
Cspec4$run(n.iter)
model4_samples <- as.matrix(Cspec4$mvSamples)

# (9) Save results
save(model4_samples, file = "model_output/model4_results_nimble.RData")
cat("Model 4 completed and saved.\n")

cat("\nAll done: Models 3 & 4 with loops for crossbasis in NIMBLE.\n")
