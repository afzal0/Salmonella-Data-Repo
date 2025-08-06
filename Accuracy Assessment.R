###############################################################################
# COMPREHENSIVE MODEL ACCURACY ASSESSMENT
# Evaluates accuracy of Models 3 & 4 using various metrics
###############################################################################

# Load libraries
library(nimble)
library(dplyr)
library(ggplot2)
library(coda)
library(reshape2)
library(lubridate)
library(gridExtra)

###############################################################################
# 1) LOAD DATA AND MODEL RESULTS
###############################################################################
cat("\n--- LOADING MODEL RESULTS AND DATA ---\n")

# Load required data
tryCatch({
  # Case-crossover data for Model 3
  cat("Loading data_casecrossover.RData...\n")
  load("model_input/data_casecrossover.RData")
  
  # Fix duplicate column names if needed
  if(anyDuplicated(names(data_cco))) {
    cat("Fixing duplicate column names in data_cco...\n")
    names(data_cco) <- make.names(names(data_cco), unique=TRUE)
  }
  
  # Time-series data for Model 4
  cat("Loading data_timeseries.RData...\n")
  load("model_input/data_timeseries.RData")
  
  # Fix duplicate column names if needed
  if(anyDuplicated(names(data_ts))) {
    cat("Fixing duplicate column names in data_ts...\n")
    names(data_ts) <- make.names(names(data_ts), unique=TRUE)
  }
  
  # Load model results
  cat("Loading model3_results_nimble.RData...\n")
  load("model_output/model3_results_nimble.RData")
  
  cat("Loading model4_results_nimble.RData...\n")
  load("model_output/model4_results_nimble.RData")
  
  # Load DIC results if available
  cat("Loading DIC results...\n")
  dic_file_found <- FALSE
  
  # Try loading the fixed DIC results first
  if(file.exists("model_output/dic_results_fixed.RData")) {
    load("model_output/dic_results_fixed.RData")
    dic_file_found <- TRUE
    cat("Loaded DIC results from model_output/dic_results_fixed.RData\n")
  } else if(file.exists("model_output/dic_results.RData")) {
    load("model_output/dic_results.RData")
    dic_file_found <- TRUE
    cat("Loaded DIC results from model_output/dic_results.RData\n")
  } else {
    cat("No DIC results file found. DIC metrics will not be available.\n")
    # Create empty dic_results object
    dic_results <- list(model3 = NULL, model4 = NULL)
  }
  
}, error = function(e) {
  cat("ERROR loading data or model results:", conditionMessage(e), "\n")
  stop("Cannot proceed without data and model results.")
})

###############################################################################
# 2) METRIC CALCULATION FUNCTIONS
###############################################################################

# Function to calculate relative metrics
calculate_relative_metrics <- function(observed, predicted) {
  # Remove zeros and NA values from observed data
  valid_idx <- which(observed > 0 & !is.na(observed) & !is.na(predicted))
  
  if(length(valid_idx) == 0) {
    return(list(MAPE = NA, SMAPE = NA))
  }
  
  observed <- observed[valid_idx]
  predicted <- predicted[valid_idx]
  
  # Calculate relative metrics
  mape <- mean(abs((observed - predicted)/observed)) * 100
  smape <- mean(200 * abs(observed - predicted)/(abs(observed) + abs(predicted)))
  
  return(list(MAPE = mape, SMAPE = smape))
}

# Function to calculate autocorrelations
calculate_autocorrelations <- function(observed, predicted) {
  # Calculate residuals
  residuals <- observed - predicted
  
  # Calculate autocorrelations at lag 1
  if(length(residuals) < 3) {
    return(list(QAC = NA, RAC = NA))
  }
  
  ac_residuals <- acf(residuals, lag.max = 1, plot = FALSE)$acf[2]
  ac_observed <- acf(observed, lag.max = 1, plot = FALSE)$acf[2]
  
  # Calculate quality of autocorrelation ratio (QAC)
  qac <- abs(ac_residuals / ac_observed)
  
  # Calculate relative autocorrelation (RAC)
  rac <- ac_residuals / ac_observed
  
  return(list(QAC = qac, RAC = rac))
}

# Function to calculate all metrics
calculate_all_metrics <- function(observed, predicted) {
  # Remove NA values
  valid_idx <- complete.cases(observed, predicted)
  observed <- observed[valid_idx]
  predicted <- predicted[valid_idx]
  
  if(length(observed) == 0) {
    return(list(
      RMSE = NA, NRMSE = NA, MAE = NA, R2 = NA, Adj_R2 = NA,
      MAPE = NA, SMAPE = NA, AIC = NA, BIC = NA, Deviance = NA,
      QAC = NA, RAC = NA, QAIC = NA
    ))
  }
  
  # Basic metrics
  rmse <- sqrt(mean((observed - predicted)^2))
  mae <- mean(abs(observed - predicted))
  r2 <- cor(observed, predicted)^2
  
  # Relative metrics with protection against zeros
  rel_metrics <- calculate_relative_metrics(observed, predicted)
  
  # Calculate normalized RMSE
  nrmse <- rmse / mean(observed)
  
  # Calculate adjusted R-squared
  n <- length(observed)
  p <- 4  # number of predictors (climate variables)
  adj_r2 <- 1 - (1 - r2) * (n - 1)/(n - p - 1)
  
  # Calculate AIC for Poisson model
  # Use pmax to ensure lambda values are positive
  pos_predicted <- pmax(predicted, 1e-10)
  deviance <- -2 * sum(dpois(observed, pos_predicted, log=TRUE))
  aic <- deviance + 2 * p
  
  # Calculate BIC
  bic <- deviance + log(n) * p
  
  # Calculate QAIC (Quasi-AIC)
  # Estimate overdispersion parameter (ĉ)
  residuals <- observed - predicted
  phi <- sum(residuals^2) / (n - p)
  c_hat <- max(1, phi)  # Use at least 1 to avoid reducing penalty
  qaic <- deviance/c_hat + 2 * p
  
  # Calculate autocorrelations
  ac_metrics <- calculate_autocorrelations(observed, predicted)
  
  # Combine all metrics
  metrics <- list(
    RMSE = rmse,
    NRMSE = nrmse,
    MAE = mae,
    R2 = r2,
    Adj_R2 = adj_r2,
    MAPE = rel_metrics$MAPE,
    SMAPE = rel_metrics$SMAPE,
    AIC = aic,
    BIC = bic,
    QAIC = qaic,
    Deviance = deviance,
    QAC = ac_metrics$QAC,
    RAC = ac_metrics$RAC,
    C_hat = c_hat  # Also return the overdispersion parameter
  )
  
  return(metrics)
}

# Function to calculate metrics by case load
calculate_load_specific_metrics <- function(observed, predicted) {
  # Create categories based on case load
  quantiles <- quantile(observed, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  # Ensure unique breakpoints
  if(length(unique(quantiles)) < 5) {
    # If not enough unique quantiles, use a simple approach
    categories <- cut(observed, 
                      breaks = c(-0.1, 5, 10, 20, max(observed, na.rm = TRUE)), 
                      labels = c("Low", "Medium-Low", "Medium-High", "High"),
                      include.lowest = TRUE)
  } else {
    categories <- cut(observed, breaks = quantiles, 
                      labels = c("Low", "Medium-Low", "Medium-High", "High"),
                      include.lowest = TRUE)
  }
  
  # Calculate metrics for each category
  load_metrics <- by(data.frame(observed = observed, predicted = predicted),
                     categories,
                     function(x) calculate_all_metrics(x$observed, x$predicted))
  
  return(load_metrics)
}

# Function to calculate metrics by season
calculate_seasonal_metrics <- function(observed, predicted, dates) {
  # Define seasons (for Southern Hemisphere, adjust as needed)
  months <- month(dates)
  seasons <- cut(months, breaks = c(0,3,6,9,12),
                 labels = c("Summer", "Autumn", "Winter", "Spring"))
  
  # Calculate metrics for each season
  seasonal_metrics <- by(data.frame(observed = observed, predicted = predicted),
                         seasons,
                         function(x) calculate_all_metrics(x$observed, x$predicted))
  
  return(seasonal_metrics)
}

# Function to create performance plots
create_performance_plots <- function(observed, predicted, dates, region = NULL) {
  # Create data frame for plotting
  plot_data <- data.frame(
    Date = dates,
    Observed = observed,
    Predicted = predicted,
    Residuals = observed - predicted
  )
  
  # Create base title
  base_title <- ifelse(is.null(region), "Overall", paste("LHD:", region))
  
  # Scatter plot
  p1 <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste(base_title, "- Predicted vs Observed"))
  
  # Time series plot
  p2 <- ggplot(plot_data) +
    geom_line(aes(x = Date, y = Observed, color = "Observed")) +
    geom_line(aes(x = Date, y = Predicted, color = "Predicted")) +
    scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
    theme_minimal() +
    labs(title = paste(base_title, "- Time Series"), color = "")
  
  # Residual plot
  p3 <- ggplot(plot_data, aes(x = Predicted, y = Residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste(base_title, "- Residuals"))
  
  # QQ plot
  p4 <- ggplot(plot_data, aes(sample = Residuals)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal() +
    labs(title = paste(base_title, "- Normal Q-Q Plot"))
  
  return(list(scatter = p1, timeseries = p2, residuals = p3, qq = p4))
}

# Function to print metrics nicely
print_metrics_summary <- function(metrics) {
  # Check for NA values
  if(is.na(metrics$RMSE)) {
    return(data.frame(
      RMSE = "NA", NRMSE = "NA", R2 = "NA", SMAPE = "NA", AIC = "NA"
    ))
  }
  
  data.frame(
    RMSE = round(metrics$RMSE, 2),
    NRMSE = round(metrics$NRMSE, 2),
    R2 = round(metrics$R2, 3),
    SMAPE = round(metrics$SMAPE, 2),
    AIC = round(metrics$AIC, 1)
  )
}

###############################################################################
# 3) EXTRACT MODEL PREDICTIONS
###############################################################################
cat("\n--- EXTRACTING MODEL PREDICTIONS ---\n")

extract_predictions <- function(model_samples, model_name) {
  # Extract lambda columns
  lambda_cols <- grep("^lambda\\[", colnames(model_samples), value=TRUE)
  
  if(length(lambda_cols) == 0) {
    cat("ERROR: No lambda parameters found in", model_name, "samples.\n")
    return(NULL)
  }
  
  # Calculate posterior mean for each lambda
  lambda_means <- colMeans(model_samples[, lambda_cols], na.rm = TRUE)
  
  # Check for NA values
  na_count <- sum(is.na(lambda_means))
  if(na_count > 0) {
    cat("WARNING:", na_count, "NA values found in predicted values for", model_name, "\n")
  }
  
  return(lambda_means)
}

# Extract predictions for Model 3 (Case-Crossover)
model3_predictions <- extract_predictions(model3_samples, "Model 3")
cat("Model 3:", length(model3_predictions), "predictions extracted\n")

# Extract predictions for Model 4 (Time-Series)
model4_predictions <- extract_predictions(model4_samples, "Model 4")
cat("Model 4:", length(model4_predictions), "predictions extracted\n")

###############################################################################
# 4) PREPARE DATA FOR ACCURACY ASSESSMENT
###############################################################################
cat("\n--- PREPARING DATA FOR ACCURACY ASSESSMENT ---\n")

# Ensure data has proper date format
if(is.character(data_cco$Date)) {
  data_cco$Date <- as.Date(data_cco$Date)
}
if(is.character(data_ts$Date)) {
  data_ts$Date <- as.Date(data_ts$Date)
}

# Create analysis data frames with observed and predicted values
model3_analysis <- data.frame(
  Date = data_cco$Date,
  LHD_code = data_cco$LHD_code,
  Observed = data_cco$Salmonella,
  Predicted = model3_predictions
)

model4_analysis <- data.frame(
  Date = data_ts$Date,
  LHD_code = data_ts$LHD_code,
  Observed = data_ts$Salmonella,
  Predicted = model4_predictions
)

###############################################################################
# 5) CALCULATE OVERALL ACCURACY METRICS
###############################################################################
cat("\n--- CALCULATING OVERALL ACCURACY METRICS ---\n")

# Calculate overall metrics for Model 3
model3_overall_metrics <- calculate_all_metrics(
  model3_analysis$Observed, 
  model3_analysis$Predicted
)

# Calculate overall metrics for Model 4
model4_overall_metrics <- calculate_all_metrics(
  model4_analysis$Observed, 
  model4_analysis$Predicted
)

# Print overall metrics
cat("\nModel 3 (Case-Crossover) Overall Metrics:\n")
print(print_metrics_summary(model3_overall_metrics))

# Add DIC metrics if available
if(!is.null(dic_results$model3)) {
  cat("\nModel 3 DIC-related Metrics:\n")
  cat("DIC:", ifelse(!is.na(dic_results$model3$DIC), round(dic_results$model3$DIC, 2), "NA"), "\n")
  cat("QAIC:", ifelse(!is.na(dic_results$model3$QAIC), round(dic_results$model3$QAIC, 2), "NA"), "\n")
  cat("Effective Parameters (pD):", ifelse(!is.na(dic_results$model3$pD), round(dic_results$model3$pD, 2), "NA"), "\n")
  cat("Deviance (Dbar):", ifelse(!is.na(dic_results$model3$Dbar), round(dic_results$model3$Dbar, 2), "NA"), "\n")
  cat("Overdispersion (from fit): ", round(model3_overall_metrics$C_hat, 3), "\n")
}

cat("\nModel 4 (Time-Series) Overall Metrics:\n")
print(print_metrics_summary(model4_overall_metrics))

# Add DIC metrics if available
if(!is.null(dic_results$model4)) {
  cat("\nModel 4 DIC-related Metrics:\n")
  cat("DIC:", ifelse(!is.na(dic_results$model4$DIC), round(dic_results$model4$DIC, 2), "NA"), "\n")
  cat("QAIC:", ifelse(!is.na(dic_results$model4$QAIC), round(dic_results$model4$QAIC, 2), "NA"), "\n")
  cat("Effective Parameters (pD):", ifelse(!is.na(dic_results$model4$pD), round(dic_results$model4$pD, 2), "NA"), "\n")
  cat("Deviance (Dbar):", ifelse(!is.na(dic_results$model4$Dbar), round(dic_results$model4$Dbar, 2), "NA"), "\n")
  cat("Overdispersion (from fit): ", round(model4_overall_metrics$C_hat, 3), "\n")
}

###############################################################################
# 6) CALCULATE REGION-SPECIFIC METRICS
###############################################################################
cat("\n--- CALCULATING REGION-SPECIFIC METRICS ---\n")

# Calculate metrics by LHD for Model 3
model3_lhd_metrics <- list()
for(lhd in unique(model3_analysis$LHD_code)) {
  lhd_data <- model3_analysis[model3_analysis$LHD_code == lhd, ]
  model3_lhd_metrics[[as.character(lhd)]] <- calculate_all_metrics(
    lhd_data$Observed, 
    lhd_data$Predicted
  )
}

# Calculate metrics by LHD for Model 4
model4_lhd_metrics <- list()
for(lhd in unique(model4_analysis$LHD_code)) {
  lhd_data <- model4_analysis[model4_analysis$LHD_code == lhd, ]
  model4_lhd_metrics[[as.character(lhd)]] <- calculate_all_metrics(
    lhd_data$Observed, 
    lhd_data$Predicted
  )
}

# Print region-specific metrics
cat("\nModel 3 (Case-Crossover) Performance by LHD:\n")
model3_lhd_summary <- do.call(rbind, lapply(names(model3_lhd_metrics), function(lhd) {
  data.frame(LHD = lhd, print_metrics_summary(model3_lhd_metrics[[lhd]]))
}))
print(model3_lhd_summary)

cat("\nModel 4 (Time-Series) Performance by LHD:\n")
model4_lhd_summary <- do.call(rbind, lapply(names(model4_lhd_metrics), function(lhd) {
  data.frame(LHD = lhd, print_metrics_summary(model4_lhd_metrics[[lhd]]))
}))
print(model4_lhd_summary)

###############################################################################
# 7) CALCULATE LOAD-SPECIFIC METRICS
###############################################################################
cat("\n--- CALCULATING LOAD-SPECIFIC METRICS ---\n")

# Calculate metrics by case load for Model 3
model3_load_metrics <- calculate_load_specific_metrics(
  model3_analysis$Observed, 
  model3_analysis$Predicted
)

# Calculate metrics by case load for Model 4
model4_load_metrics <- calculate_load_specific_metrics(
  model4_analysis$Observed, 
  model4_analysis$Predicted
)

# Print load-specific metrics
cat("\nModel 3 (Case-Crossover) Performance by Case Load:\n")
for(level in names(model3_load_metrics)) {
  cat("\n  Level:", level, "\n")
  print(print_metrics_summary(model3_load_metrics[[level]]))
}

cat("\nModel 4 (Time-Series) Performance by Case Load:\n")
for(level in names(model4_load_metrics)) {
  cat("\n  Level:", level, "\n")
  print(print_metrics_summary(model4_load_metrics[[level]]))
}

###############################################################################
# 8) CALCULATE SEASONAL METRICS
###############################################################################
cat("\n--- CALCULATING SEASONAL METRICS ---\n")

# Calculate metrics by season for Model 3
model3_seasonal_metrics <- calculate_seasonal_metrics(
  model3_analysis$Observed, 
  model3_analysis$Predicted,
  model3_analysis$Date
)

# Calculate metrics by season for Model 4
model4_seasonal_metrics <- calculate_seasonal_metrics(
  model4_analysis$Observed, 
  model4_analysis$Predicted,
  model4_analysis$Date
)

# Print seasonal metrics
cat("\nModel 3 (Case-Crossover) Performance by Season:\n")
for(season in names(model3_seasonal_metrics)) {
  cat("\n  Season:", season, "\n")
  print(print_metrics_summary(model3_seasonal_metrics[[season]]))
}

cat("\nModel 4 (Time-Series) Performance by Season:\n")
for(season in names(model4_seasonal_metrics)) {
  cat("\n  Season:", season, "\n")
  print(print_metrics_summary(model4_seasonal_metrics[[season]]))
}

###############################################################################
# 9) CREATE VISUALIZATIONS
###############################################################################
cat("\n--- CREATING VISUALIZATIONS ---\n")

# Create output directory if it doesn't exist
if(!dir.exists("model_output/accuracy_plots")) {
  dir.create("model_output/accuracy_plots", recursive = TRUE)
}

# Create overall performance plots
model3_overall_plots <- create_performance_plots(
  model3_analysis$Observed,
  model3_analysis$Predicted,
  model3_analysis$Date
)

model4_overall_plots <- create_performance_plots(
  model4_analysis$Observed,
  model4_analysis$Predicted,
  model4_analysis$Date
)

# Save overall plots
pdf("model_output/accuracy_plots/model3_overall_performance.pdf", width = 12, height = 10)
grid.arrange(
  model3_overall_plots$scatter,
  model3_overall_plots$timeseries,
  model3_overall_plots$residuals,
  model3_overall_plots$qq,
  nrow = 2
)
dev.off()

pdf("model_output/accuracy_plots/model4_overall_performance.pdf", width = 12, height = 10)
grid.arrange(
  model4_overall_plots$scatter,
  model4_overall_plots$timeseries,
  model4_overall_plots$residuals,
  model4_overall_plots$qq,
  nrow = 2
)
dev.off()

# Create region-specific plots (just for first few regions to avoid too many files)
sample_regions <- unique(model3_analysis$LHD_code)[1:min(5, length(unique(model3_analysis$LHD_code)))]

for(lhd in sample_regions) {
  # Model 3 region plots
  lhd_data3 <- model3_analysis[model3_analysis$LHD_code == lhd, ]
  lhd_plots3 <- create_performance_plots(
    lhd_data3$Observed,
    lhd_data3$Predicted,
    lhd_data3$Date,
    lhd
  )
  
  pdf(paste0("model_output/accuracy_plots/model3_lhd_", lhd, "_performance.pdf"), width = 12, height = 10)
  grid.arrange(
    lhd_plots3$scatter,
    lhd_plots3$timeseries,
    lhd_plots3$residuals,
    lhd_plots3$qq,
    nrow = 2
  )
  dev.off()
  
  # Model 4 region plots
  lhd_data4 <- model4_analysis[model4_analysis$LHD_code == lhd, ]
  lhd_plots4 <- create_performance_plots(
    lhd_data4$Observed,
    lhd_data4$Predicted,
    lhd_data4$Date,
    lhd
  )
  
  pdf(paste0("model_output/accuracy_plots/model4_lhd_", lhd, "_performance.pdf"), width = 12, height = 10)
  grid.arrange(
    lhd_plots4$scatter,
    lhd_plots4$timeseries,
    lhd_plots4$residuals,
    lhd_plots4$qq,
    nrow = 2
  )
  dev.off()
}

###############################################################################
# 10) CREATE ADDITIONAL VISUALIZATIONS
###############################################################################
cat("\n--- CREATING ADDITIONAL VISUALIZATIONS ---\n")

# Create comparison bar chart with standard metrics
comparison_plot <- ggplot(model_comparison, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Model Performance Comparison",
    subtitle = "Lower values are better for all metrics except R²",
    caption = "Note: AIC, BIC, DIC and QAIC values are divided by 1000 for better visualization"
  )

# Save comparison plot
pdf("model_output/accuracy_plots/model_comparison.pdf", width = 10, height = 8)
print(comparison_plot)
dev.off()

# 1. Create Radar/Spider Chart for model comparison
# First prepare data in the right format for radar chart
if(require(fmsb)) {
  # Extract only the metrics we want to show in radar plot
  radar_metrics <- c("RMSE", "NRMSE", "SMAPE", "1-R²") # We use 1-R² so lower is better for all
  
  # Create data frame for radar chart
  radar_data <- data.frame(
    Model3 = c(
      model3_overall_metrics$RMSE,
      model3_overall_metrics$NRMSE,
      model3_overall_metrics$SMAPE,
      1 - model3_overall_metrics$R2
    ),
    Model4 = c(
      model4_overall_metrics$RMSE,
      model4_overall_metrics$NRMSE,
      model4_overall_metrics$SMAPE,
      1 - model4_overall_metrics$R2
    )
  )
  
  rownames(radar_data) <- radar_metrics
  
  # Add max and min values required by fmsb
  radar_data <- rbind(
    apply(radar_data, 2, max) * 1.1, # Max value (slightly higher for better visualization)
    rep(0, ncol(radar_data)),        # Min value
    radar_data
  )
  
  # Create radar plot
  pdf("model_output/accuracy_plots/radar_comparison.pdf", width = 10, height = 8)
  par(mar = c(1, 1, 3, 1))
  radarchart(
    radar_data,
    pcol = c("red", "blue"),
    pfcol = c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.3)),
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    title = "Model Comparison (Lower is Better)"
  )
  legend(
    "topright",
    legend = c("Model 3 (Case-Crossover)", "Model 4 (Time-Series)"),
    col = c("red", "blue"),
    lty = 1,
    lwd = 2,
    pch = 16,
    bty = "n"
  )
  dev.off()
}

# 2. Create heatmap of performance by LHD
# Extract RMSE for each LHD
lhd_rmse_data <- data.frame(
  LHD = numeric(),
  Model = character(),
  RMSE = numeric(),
  NRMSE = numeric(),
  R2 = numeric()
)

# Collect metrics by LHD
for(lhd in names(model3_lhd_metrics)) {
  lhd_rmse_data <- rbind(lhd_rmse_data, data.frame(
    LHD = as.numeric(lhd),
    Model = "Model 3",
    RMSE = model3_lhd_metrics[[lhd]]$RMSE,
    NRMSE = model3_lhd_metrics[[lhd]]$NRMSE,
    R2 = model3_lhd_metrics[[lhd]]$R2
  ))
}

for(lhd in names(model4_lhd_metrics)) {
  lhd_rmse_data <- rbind(lhd_rmse_data, data.frame(
    LHD = as.numeric(lhd),
    Model = "Model 4",
    RMSE = model4_lhd_metrics[[lhd]]$RMSE,
    NRMSE = model4_lhd_metrics[[lhd]]$NRMSE,
    R2 = model4_lhd_metrics[[lhd]]$R2
  ))
}

# Create heatmap plot for RMSE by LHD
heatmap_rmse <- ggplot(lhd_rmse_data, aes(x = factor(LHD), y = Model, fill = RMSE)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(lhd_rmse_data$RMSE, na.rm = TRUE)) +
  theme_minimal() +
  labs(
    title = "RMSE by Local Health District",
    x = "Local Health District (LHD) Code",
    y = "Model",
    fill = "RMSE"
  ) +
  theme(axis.text.x = element_text(angle = 0))

# Create heatmap plot for R² by LHD
heatmap_r2 <- ggplot(lhd_rmse_data, aes(x = factor(LHD), y = Model, fill = R2)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(lhd_rmse_data$R2, na.rm = TRUE)) +
  theme_minimal() +
  labs(
    title = "R² by Local Health District",
    x = "Local Health District (LHD) Code",
    y = "Model",
    fill = "R²"
  ) +
  theme(axis.text.x = element_text(angle = 0))

# Save heatmap plots
pdf("model_output/accuracy_plots/lhd_heatmap.pdf", width = 12, height = 6)
grid.arrange(heatmap_rmse, heatmap_r2, ncol = 1)
dev.off()

# 3. Create seasonal performance visualization
# Extract seasonal metrics
seasonal_data <- data.frame(
  Season = character(),
  Model = character(),
  RMSE = numeric(),
  R2 = numeric()
)

# Process Model 3 seasonal metrics
for(season in names(model3_seasonal_metrics)) {
  seasonal_data <- rbind(seasonal_data, data.frame(
    Season = season,
    Model = "Model 3",
    RMSE = model3_seasonal_metrics[[season]]$RMSE,
    R2 = model3_seasonal_metrics[[season]]$R2
  ))
}

# Process Model 4 seasonal metrics
for(season in names(model4_seasonal_metrics)) {
  seasonal_data <- rbind(seasonal_data, data.frame(
    Season = season,
    Model = "Model 4",
    RMSE = model4_seasonal_metrics[[season]]$RMSE,
    R2 = model4_seasonal_metrics[[season]]$R2
  ))
}

# Order seasons correctly
seasonal_data$Season <- factor(seasonal_data$Season, 
                               levels = c("Summer", "Autumn", "Winter", "Spring"))

# Create seasonal performance bar plots
seasonal_rmse <- ggplot(seasonal_data, aes(x = Season, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Seasonal RMSE Comparison",
    subtitle = "Lower values indicate better performance",
    x = "Season",
    y = "RMSE"
  )

seasonal_r2 <- ggplot(seasonal_data, aes(x = Season, y = R2, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Seasonal R² Comparison",
    subtitle = "Higher values indicate better performance",
    x = "Season",
    y = "R²"
  )

# Save seasonal plots
pdf("model_output/accuracy_plots/seasonal_performance.pdf", width = 12, height = 8)
grid.arrange(seasonal_rmse, seasonal_r2, ncol = 1)
dev.off()

# 4. Create case load performance visualization
# Extract case load metrics
load_data <- data.frame(
  Level = character(),
  Model = character(),
  RMSE = numeric(),
  R2 = numeric()
)

# Process Model 3 load metrics
for(level in names(model3_load_metrics)) {
  load_data <- rbind(load_data, data.frame(
    Level = level,
    Model = "Model 3",
    RMSE = model3_load_metrics[[level]]$RMSE,
    R2 = model3_load_metrics[[level]]$R2
  ))
}

# Process Model 4 load metrics
for(level in names(model4_load_metrics)) {
  load_data <- rbind(load_data, data.frame(
    Level = level,
    Model = "Model 4",
    RMSE = model4_load_metrics[[level]]$RMSE,
    R2 = model4_load_metrics[[level]]$R2
  ))
}

# Order load levels correctly
load_data$Level <- factor(load_data$Level, 
                          levels = c("Low", "Medium-Low", "Medium-High", "High"))

# Create case load performance bar plots
load_rmse <- ggplot(load_data, aes(x = Level, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Case Load RMSE Comparison",
    subtitle = "Lower values indicate better performance",
    x = "Case Load Level",
    y = "RMSE"
  )

load_r2 <- ggplot(load_data, aes(x = Level, y = R2, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Case Load R² Comparison",
    subtitle = "Higher values indicate better performance",
    x = "Case Load Level",
    y = "R²"
  )

# Save case load plots
pdf("model_output/accuracy_plots/caseload_performance.pdf", width = 12, height = 8)
grid.arrange(load_rmse, load_r2, ncol = 1)
dev.off()

# 5. Create information criteria comparison
# Extract information criteria from both models
ic_data <- data.frame(
  Metric = c("AIC", "BIC", "QAIC", "DIC (MCMC)", "QAIC (MCMC)"),
  Model3 = c(
    model3_overall_metrics$AIC,
    model3_overall_metrics$BIC,
    model3_overall_metrics$QAIC,
    ifelse(!is.null(dic_results$model3) && !is.na(dic_results$model3$DIC), 
           dic_results$model3$DIC, NA),
    ifelse(!is.null(dic_results$model3) && !is.na(dic_results$model3$QAIC), 
           dic_results$model3$QAIC, NA)
  ),
  Model4 = c(
    model4_overall_metrics$AIC,
    model4_overall_metrics$BIC,
    model4_overall_metrics$QAIC,
    ifelse(!is.null(dic_results$model4) && !is.na(dic_results$model4$DIC), 
           dic_results$model4$DIC, NA),
    ifelse(!is.null(dic_results$model4) && !is.na(dic_results$model4$QAIC), 
           dic_results$model4$QAIC, NA)
  )
)

# Convert to long format for plotting
ic_data_long <- reshape2::melt(ic_data, id.vars = "Metric", 
                               variable.name = "Model", value.name = "Value")

# Create information criteria bar plot
ic_plot <- ggplot(ic_data_long, aes(x = Metric, y = Value/1000, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Information Criteria Comparison",
    subtitle = "Lower values indicate better model fit",
    x = "Information Criterion",
    y = "Value (thousands)",
    caption = "DIC and QAIC (MCMC) are calculated from MCMC samples"
  )

# Save information criteria plot
pdf("model_output/accuracy_plots/information_criteria.pdf", width = 10, height = 8)
print(ic_plot)
dev.off()

# 6. Create residual patterns by month
# Extract month from dates
model3_analysis$Month <- month(model3_analysis$Date, label = TRUE)
model4_analysis$Month <- month(model4_analysis$Date, label = TRUE)

# Calculate residuals
model3_analysis$Residual <- model3_analysis$Observed - model3_analysis$Predicted
model4_analysis$Residual <- model4_analysis$Observed - model4_analysis$Predicted

# Create residual boxplots by month
month_residuals3 <- ggplot(model3_analysis, aes(x = Month, y = Residual)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Model 3: Residuals by Month",
    x = "Month",
    y = "Residual (Observed - Predicted)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

month_residuals4 <- ggplot(model4_analysis, aes(x = Month, y = Residual)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(
    title = "Model 4: Residuals by Month",
    x = "Month",
    y = "Residual (Observed - Predicted)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

# Save monthly residual plots
pdf("model_output/accuracy_plots/monthly_residuals.pdf", width = 12, height = 8)
grid.arrange(month_residuals3, month_residuals4, ncol = 1)
dev.off()

###############################################################################
# 11) CREATE MODEL COMPARISON PLOT
###############################################################################
cat("\n--- CREATING MODEL COMPARISON VISUALIZATIONS ---\n")

# Create a data frame for model comparison
# Include standard metrics
metrics_to_include <- c("RMSE", "NRMSE", "R²", "SMAPE", "AIC", "BIC", "QAIC")
metrics_values_model3 <- c(
  model3_overall_metrics$RMSE,
  model3_overall_metrics$NRMSE,
  model3_overall_metrics$R2,
  model3_overall_metrics$SMAPE,
  model3_overall_metrics$AIC/1000, # Scaled for better visualization
  model3_overall_metrics$BIC/1000, # Scaled for better visualization
  model3_overall_metrics$QAIC/1000 # Scaled for better visualization
)
metrics_values_model4 <- c(
  model4_overall_metrics$RMSE,
  model4_overall_metrics$NRMSE,
  model4_overall_metrics$R2,
  model4_overall_metrics$SMAPE,
  model4_overall_metrics$AIC/1000, # Scaled for better visualization
  model4_overall_metrics$BIC/1000, # Scaled for better visualization
  model4_overall_metrics$QAIC/1000 # Scaled for better visualization
)

# Add DIC metrics if available
if(exists("dic_results") && !is.null(dic_results$model3) && !is.null(dic_results$model4) &&
   !is.na(dic_results$model3$DIC) && !is.na(dic_results$model4$DIC)) {
  metrics_to_include <- c(metrics_to_include, "DIC (MCMC)")
  metrics_values_model3 <- c(metrics_values_model3, dic_results$model3$DIC/1000) # Scaled
  metrics_values_model4 <- c(metrics_values_model4, dic_results$model4$DIC/1000) # Scaled
  
  # Also add MCMC QAIC if available
  if(!is.na(dic_results$model3$QAIC) && !is.na(dic_results$model4$QAIC)) {
    metrics_to_include <- c(metrics_to_include, "QAIC (MCMC)")
    metrics_values_model3 <- c(metrics_values_model3, dic_results$model3$QAIC/1000) # Scaled
    metrics_values_model4 <- c(metrics_values_model4, dic_results$model4$QAIC/1000) # Scaled
  }
}

# Combine into a data frame
model_comparison <- data.frame(
  Model = rep(c("Model 3 (Case-Crossover)", "Model 4 (Time-Series)"), each = length(metrics_to_include)),
  Metric = rep(metrics_to_include, 2),
  Value = c(metrics_values_model3, metrics_values_model4)
)

# Create comparison plot
comparison_plot <- ggplot(model_comparison, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Model Performance Comparison",
    subtitle = "Lower values are better for all metrics except R²",
    caption = "Note: AIC and BIC values are divided by 1000 for better visualization"
  )

# Save comparison plot
pdf("model_output/accuracy_plots/model_comparison.pdf", width = 10, height = 8)
print(comparison_plot)
dev.off()

###############################################################################
# 11) SAVE RESULTS
###############################################################################
cat("\n--- SAVING ACCURACY ASSESSMENT RESULTS ---\n")

accuracy_results <- list(
  model3 = list(
    overall = model3_overall_metrics,
    by_region = model3_lhd_metrics,
    by_load = model3_load_metrics,
    by_season = model3_seasonal_metrics,
    predictions = model3_analysis
  ),
  model4 = list(
    overall = model4_overall_metrics,
    by_region = model4_lhd_metrics,
    by_load = model4_load_metrics,
    by_season = model4_seasonal_metrics,
    predictions = model4_analysis
  )
)

save(accuracy_results, file = "model_output/accuracy_assessment_results.RData")
cat("Results saved to model_output/accuracy_assessment_results.RData\n")

###############################################################################
# 12) GENERATE SUMMARY TABLE
###############################################################################
cat("\n--- GENERATING SUMMARY TABLE ---\n")

summary_table <- data.frame(
  Metric = c("RMSE", "NRMSE", "R²", "Adj. R²", "MAE", "SMAPE", "AIC", "BIC", "QAIC"),
  Model3 = c(
    round(model3_overall_metrics$RMSE, 2),
    round(model3_overall_metrics$NRMSE, 3),
    round(model3_overall_metrics$R2, 3),
    round(model3_overall_metrics$Adj_R2, 3),
    round(model3_overall_metrics$MAE, 2),
    round(model3_overall_metrics$SMAPE, 2),
    round(model3_overall_metrics$AIC, 2),
    round(model3_overall_metrics$BIC, 2),
    round(model3_overall_metrics$QAIC, 2)
  ),
  Model4 = c(
    round(model4_overall_metrics$RMSE, 2),
    round(model4_overall_metrics$NRMSE, 3),
    round(model4_overall_metrics$R2, 3),
    round(model4_overall_metrics$Adj_R2, 3),
    round(model4_overall_metrics$MAE, 2),
    round(model4_overall_metrics$SMAPE, 2),
    round(model4_overall_metrics$AIC, 2),
    round(model4_overall_metrics$BIC, 2),
    round(model4_overall_metrics$QAIC, 2)
  ),
  BestModel = c(
    ifelse(model3_overall_metrics$RMSE < model4_overall_metrics$RMSE, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$NRMSE < model4_overall_metrics$NRMSE, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$R2 > model4_overall_metrics$R2, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$Adj_R2 > model4_overall_metrics$Adj_R2, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$MAE < model4_overall_metrics$MAE, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$SMAPE < model4_overall_metrics$SMAPE, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$AIC < model4_overall_metrics$AIC, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$BIC < model4_overall_metrics$BIC, "Model 3", "Model 4"),
    ifelse(model3_overall_metrics$QAIC < model4_overall_metrics$QAIC, "Model 3", "Model 4")
  )
)

# Update summary table with DIC/QAIC values
if(exists("dic_results") && !is.null(dic_results$model3) && !is.null(dic_results$model4)) {
  # Append DIC metrics
  dic_rows <- data.frame(
    Metric = c("DIC (MCMC)", "QAIC (MCMC)", "QAIC (Fit)", "Overdispersion"),
    Model3 = c(
      ifelse(!is.na(dic_results$model3$DIC), round(dic_results$model3$DIC, 2), "NA"),
      ifelse(!is.na(dic_results$model3$QAIC), round(dic_results$model3$QAIC, 2), "NA"),
      round(model3_overall_metrics$QAIC, 2),
      round(model3_overall_metrics$C_hat, 3)
    ),
    Model4 = c(
      ifelse(!is.na(dic_results$model4$DIC), round(dic_results$model4$DIC, 2), "NA"),
      ifelse(!is.na(dic_results$model4$QAIC), round(dic_results$model4$QAIC, 2), "NA"),
      round(model4_overall_metrics$QAIC, 2),
      round(model4_overall_metrics$C_hat, 3)
    ),
    BestModel = c(
      ifelse(!is.na(dic_results$model3$DIC) && !is.na(dic_results$model4$DIC),
             ifelse(dic_results$model3$DIC < dic_results$model4$DIC, "Model 3", "Model 4"), "NA"),
      ifelse(!is.na(dic_results$model3$QAIC) && !is.na(dic_results$model4$QAIC),
             ifelse(dic_results$model3$QAIC < dic_results$model4$QAIC, "Model 3", "Model 4"), "NA"),
      ifelse(model3_overall_metrics$QAIC < model4_overall_metrics$QAIC, "Model 3", "Model 4"),
      "N/A"
    )
  )
  
  # Combine with existing summary table
  summary_table <- rbind(summary_table, dic_rows)
}

# Print and save summary table
print(summary_table)
write.csv(summary_table, "model_output/model_accuracy_summary.csv", row.names = FALSE)

# Generate comprehensive model comparison table as HTML
if(require(knitr) && require(kableExtra)) {
  cat("\nGenerating HTML version of model comparison table...\n")
  
  # Style the table with colors to highlight best model
  styled_table <- kable(summary_table, format = "html", caption = "Model Comparison Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    column_spec(4, background = ifelse(summary_table$BestModel == "Model 3", "#d4edda", 
                                       ifelse(summary_table$BestModel == "Model 4", "#f8d7da", "white")))
  
  # Save to file
  writeLines(styled_table, "model_output/model_comparison_table.html")
  cat("Saved HTML table to model_output/model_comparison_table.html\n")
}

cat("\nAccuracy assessment complete. Results saved to model_output directory.\n")