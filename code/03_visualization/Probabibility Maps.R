###############################################################################
# relative_risk_probability_maps.R
#
# Creates maps showing probability of exceeding relative risk thresholds:
#  - RR > 1.0
#  - RR > 1.2
#  - RR > 1.5
#
# For Model 3 and Model 4 and all variables (Temperature, Rainfall, Flooding, Humidity)
# in a 3x3 grid format with consistent color scheme
#
###############################################################################

#-----------------------------
# 1) Load libraries
#-----------------------------
library(sf)          
library(sp)
library(spdep)       
library(ggplot2)     
library(dplyr)       
library(coda)        
library(patchwork)   
library(scales)      # For number formatting
library(gridExtra)   # For grid.arrange
library(RColorBrewer) # For color palettes

#-----------------------------
# 2) Load shapefile & MCMC results
#-----------------------------
# Load shapefile (adjust path if needed)
shp_path <- "model_input/NSW_LHD_Boundaries.shp"
if(!file.exists(shp_path)) {
  shp_path <- file.path("data/raw/NSW_LHD_Boundaries.shp")
  if(!file.exists(shp_path)) {
    stop("Shapefile not found. Please provide the correct path.")
  }
}

shp_sf <- st_read(shp_path)
shp_sf$LHD_code <- as.integer(factor(shp_sf$lhd_name))
lhd_names <- shp_sf$lhd_name

# Spatial adjacency list
shp_sp <- as(shp_sf, "Spatial")
nb_list <- poly2nb(shp_sp, row.names=shp_sp@data$LHD_code)

# Load NIMBLE model results
cat("Loading NIMBLE model results...\n")
load("model_output/model3_results_nimble.RData")  # => model3_samples
load("model_output/model4_results_nimble.RData")  # => model4_samples

# Convert to matrix format if not already
if(!is.matrix(model3_samples)) {
  cat("Converting model3_samples to matrix...\n")
  m3_mat <- as.matrix(model3_samples)
} else {
  m3_mat <- model3_samples
}

if(!is.matrix(model4_samples)) {
  cat("Converting model4_samples to matrix...\n")
  m4_mat <- as.matrix(model4_samples)
} else {
  m4_mat <- model4_samples
}

#-----------------------------
# 3) Helper functions
#-----------------------------
# Function to extract MCMC samples for each region
getNetEffect <- function(mcmc_matrix, param_prefix, n_region) {
  region_list <- vector("list", n_region)
  for(r in seq_len(n_region)) {
    pattern_r <- paste0("^", param_prefix, "\\[", r, ",")
    cols_r    <- grep(pattern_r, colnames(mcmc_matrix), value=TRUE)
    
    if(length(cols_r) == 0) {
      # Try alternative patterns for different parameter naming conventions
      alt_patterns <- c(
        paste0("^", param_prefix, "_", r, "_"),  # For beta_temp_1_1 format
        paste0("^", param_prefix, r, "_"),       # For betatemp1_1 format
        paste0("^", param_prefix, "_", r, ",")   # For beta_temp_1,1 format
      )
      
      for(alt_pattern in alt_patterns) {
        cols_r <- grep(alt_pattern, colnames(mcmc_matrix), value=TRUE)
        if(length(cols_r) > 0) break
      }
      
      if(length(cols_r) == 0) {
        warning("No columns found for pattern: ", pattern_r, " or alternatives")
        region_list[[r]] <- rep(NA, nrow(mcmc_matrix))
        next
      }
    }
    
    # Return the sum of effects across all parameters for this region
    region_list[[r]] <- rowSums(mcmc_matrix[, cols_r, drop=FALSE])
  }
  do.call(cbind, region_list) # [draws x region]
}

# Calculate exceedance probabilities for specified threshold
calculate_exceedance_prob <- function(net_effect_matrix, threshold) {
  if(all(is.na(net_effect_matrix))) {
    return(rep(NA, ncol(net_effect_matrix)))
  }
  
  # Convert to relative risk and calculate probability of exceeding threshold
  rr_matrix <- exp(net_effect_matrix)
  exceed_probs <- colMeans(rr_matrix > threshold, na.rm = TRUE)
  
  return(exceed_probs)
}

# Function to check valid data
check_valid_data <- function(net_effect, variable_name) {
  valid_pct <- sum(!is.na(net_effect)) / length(net_effect) * 100
  cat(sprintf("%s has %.1f%% valid data\n", variable_name, valid_pct))
  return(valid_pct > 0)  # Return TRUE if we have at least some valid data
}

# Map plotting function for probability maps
plot_probability_map <- function(
    sf_object, prob_column,
    title_text = NULL,
    var_name = NULL,
    threshold = NULL
) {
  plot_title <- if(is.null(title_text)) {
    sprintf("%s: Probability of Exceeding RR=%.1f", var_name, threshold)
  } else {
    title_text
  }
  
  p <- ggplot() +
    geom_sf(
      data = sf_object,
      aes_string(fill = prob_column),
      color = "#666666",
      size = 0.1
    ) +
    scale_fill_gradientn(
      name = sprintf("Prob(RR > %.1f)", threshold),
      colors = brewer.pal(9, "YlOrRd"),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey80"
    ) +
    labs(title = plot_title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      plot.title = element_text(size = 11, face = "plain"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.position = "right",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  return(p)
}

#-----------------------------
# 4) Process MCMC samples
#-----------------------------
n_regions <- length(unique(shp_sf$LHD_code))

# MODEL 3
cat("Calculating Model 3 effect summaries...\n")
m3_temp_net  <- getNetEffect(m3_mat, "beta_temp",  n_regions)
m3_rain_net  <- getNetEffect(m3_mat, "beta_rain",  n_regions)
m3_flood_net <- getNetEffect(m3_mat, "beta_flood", n_regions)
m3_humid_net <- getNetEffect(m3_mat, "beta_humidity", n_regions)

# Check for valid data
valid_m3_temp  <- check_valid_data(m3_temp_net, "Model 3 temperature")
valid_m3_rain  <- check_valid_data(m3_rain_net, "Model 3 rainfall")
valid_m3_flood <- check_valid_data(m3_flood_net, "Model 3 flood")
valid_m3_humid <- check_valid_data(m3_humid_net, "Model 3 humidity")

# MODEL 4
cat("Calculating Model 4 effect summaries...\n")
m4_temp_net  <- getNetEffect(m4_mat, "beta_temp",  n_regions)
m4_rain_net  <- getNetEffect(m4_mat, "beta_rain",  n_regions)
m4_flood_net <- getNetEffect(m4_mat, "beta_flood", n_regions)
m4_humid_net <- getNetEffect(m4_mat, "beta_humidity", n_regions)

# Check for valid data
valid_m4_temp  <- check_valid_data(m4_temp_net, "Model 4 temperature")
valid_m4_rain  <- check_valid_data(m4_rain_net, "Model 4 rainfall")
valid_m4_flood <- check_valid_data(m4_flood_net, "Model 4 flood")
valid_m4_humid <- check_valid_data(m4_humid_net, "Model 4 humidity")

#-----------------------------
# 5) Calculate probability maps for MODEL 3
#-----------------------------
# Define thresholds
thresholds <- c(1.0, 1.2, 1.5)

# For each variable and threshold, calculate exceedance probabilities and prepare maps
m3_sf <- shp_sf

# Temperature probabilities
if(valid_m3_temp) {
  for(t in thresholds) {
    col_name <- paste0("prob_temp_", gsub("\\.", "_", t))
    m3_sf[[col_name]] <- calculate_exceedance_prob(m3_temp_net, t)
  }
}

# Rainfall probabilities
if(valid_m3_rain) {
  for(t in thresholds) {
    col_name <- paste0("prob_rain_", gsub("\\.", "_", t))
    m3_sf[[col_name]] <- calculate_exceedance_prob(m3_rain_net, t)
  }
}

# Flooding probabilities 
if(valid_m3_flood) {
  for(t in thresholds) {
    col_name <- paste0("prob_flood_", gsub("\\.", "_", t))
    m3_sf[[col_name]] <- calculate_exceedance_prob(m3_flood_net, t)
  }
}

# Humidity probabilities
if(valid_m3_humid) {
  for(t in thresholds) {
    col_name <- paste0("prob_humid_", gsub("\\.", "_", t))
    m3_sf[[col_name]] <- calculate_exceedance_prob(m3_humid_net, t)
  }
}

#-----------------------------
# 6) Calculate probability maps for MODEL 4
#-----------------------------
m4_sf <- shp_sf

# Temperature probabilities
if(valid_m4_temp) {
  for(t in thresholds) {
    col_name <- paste0("prob_temp_", gsub("\\.", "_", t))
    m4_sf[[col_name]] <- calculate_exceedance_prob(m4_temp_net, t)
  }
}

# Rainfall probabilities
if(valid_m4_rain) {
  for(t in thresholds) {
    col_name <- paste0("prob_rain_", gsub("\\.", "_", t))
    m4_sf[[col_name]] <- calculate_exceedance_prob(m4_rain_net, t)
  }
}

# Flooding probabilities 
if(valid_m4_flood) {
  for(t in thresholds) {
    col_name <- paste0("prob_flood_", gsub("\\.", "_", t))
    m4_sf[[col_name]] <- calculate_exceedance_prob(m4_flood_net, t)
  }
}

# Humidity probabilities
if(valid_m4_humid) {
  for(t in thresholds) {
    col_name <- paste0("prob_humid_", gsub("\\.", "_", t))
    m4_sf[[col_name]] <- calculate_exceedance_prob(m4_humid_net, t)
  }
}

#-----------------------------
# 7) Create probability plot grids
#-----------------------------
# Function to create a probability map grid for a model
create_model_probability_grid <- function(model_sf, model_name, 
                                          variables, var_names, thresholds,
                                          valid_flags) {
  plots <- list()
  
  for(i in seq_along(variables)) {
    var <- variables[i]
    var_name <- var_names[i]
    
    # Skip this variable if data is not valid
    if(!valid_flags[i]) {
      for(t in thresholds) {
        # Create empty plot as placeholder
        plots[[paste0(var, "_", t)]] <- ggplot() + 
          annotate("text", x = 0, y = 0, label = "No data available") +
          theme_void() +
          labs(title = sprintf("%s: Probability of Exceeding RR=%.1f", var_name, t))
      }
      next
    }
    
    # Create plots for each threshold
    for(t in thresholds) {
      col_name <- paste0("prob_", var, "_", gsub("\\.", "_", t))
      
      plots[[paste0(var, "_", t)]] <- plot_probability_map(
        model_sf,
        col_name,
        var_name = var_name,
        threshold = t
      )
    }
  }
  
  # Arrange all plots in a grid
  plots_to_arrange <- list()
  # Add plots in the correct order for grid
  for(var_idx in seq_along(variables)) {
    var <- variables[var_idx]
    for(t in thresholds) {
      key <- paste0(var, "_", t)
      if(key %in% names(plots)) {
        plots_to_arrange[[length(plots_to_arrange) + 1]] <- plots[[key]]
      }
    }
  }
  
  # Create the grid
  if(length(plots_to_arrange) > 0) {
    grid_plot <- gridExtra::grid.arrange(
      grobs = plots_to_arrange,
      ncol = length(thresholds),
      top = grid::textGrob(
        sprintf("%s: Probability Analysis", model_name),
        gp = grid::gpar(fontsize = 16, fontface = "bold"),
        just = "center"
      )
    )
    return(grid_plot)
  } else {
    return(NULL)
  }
}

#-----------------------------
# 8) Output to PDF
#-----------------------------
pdf("model_output/relative_risk_probability_maps.pdf", width = 11, height = 8.5)

# Create Model 3 probability grid
vars <- c("temp", "rain", "flood", "humid")
var_names <- c("Temperature", "Rainfall", "Flooding", "Humidity")
valid_m3_flags <- c(valid_m3_temp, valid_m3_rain, valid_m3_flood, valid_m3_humid)

cat("Creating Model 3 probability grid...\n")
model3_grid <- create_model_probability_grid(
  m3_sf, "Model3", 
  vars, var_names, thresholds,
  valid_m3_flags
)
print(model3_grid)

# Create Model 4 probability grid
valid_m4_flags <- c(valid_m4_temp, valid_m4_rain, valid_m4_flood, valid_m4_humid)

cat("Creating Model 4 probability grid...\n")
model4_grid <- create_model_probability_grid(
  m4_sf, "Model4", 
  vars, var_names, thresholds,
  valid_m4_flags
)
print(model4_grid)

dev.off()

#-----------------------------
# 9) Save individual probability maps
#-----------------------------
# Create directory for probability maps
dir.create("model_output/probability_maps", showWarnings = FALSE)

# Create a list of all individual probability maps
m3_prob_plots <- list()
m4_prob_plots <- list()

# Generate individual probability maps for Model 3
for(i in seq_along(vars)) {
  var <- vars[i]
  var_name <- var_names[i]
  
  if(valid_m3_flags[i]) {
    for(t in thresholds) {
      col_name <- paste0("prob_", var, "_", gsub("\\.", "_", t))
      plot_key <- paste0(var, "_", t)
      
      m3_prob_plots[[plot_key]] <- plot_probability_map(
        m3_sf,
        col_name,
        var_name = var_name,
        threshold = t
      )
      
      # Save individual plot
      filename <- sprintf("model_output/probability_maps/model3_%s_rr%.1f.png", var, t)
      ggsave(filename, m3_prob_plots[[plot_key]], width = 7, height = 6, dpi = 300)
    }
  }
}

# Generate individual probability maps for Model 4
for(i in seq_along(vars)) {
  var <- vars[i]
  var_name <- var_names[i]
  
  if(valid_m4_flags[i]) {
    for(t in thresholds) {
      col_name <- paste0("prob_", var, "_", gsub("\\.", "_", t))
      plot_key <- paste0(var, "_", t)
      
      m4_prob_plots[[plot_key]] <- plot_probability_map(
        m4_sf,
        col_name,
        var_name = var_name,
        threshold = t
      )
      
      # Save individual plot
      filename <- sprintf("model_output/probability_maps/model4_%s_rr%.1f.png", var, t)
      ggsave(filename, m4_prob_plots[[plot_key]], width = 7, height = 6, dpi = 300)
    }
  }
}

cat("\nProbability maps completed. Output files:\n")
cat("- PDF report: model_output/relative_risk_probability_maps.pdf\n")
cat("- Individual maps: model_output/probability_maps/\n")