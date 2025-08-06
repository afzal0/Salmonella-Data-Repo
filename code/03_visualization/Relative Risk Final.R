###############################################################################
# relative_risk_maps_fixed.R - Adjusted for NIMBLE Models
#
# Creates:
#  (1) Model 3 RR maps (Temp, Rain, Flood, Humidity) with improved readable legends
#  (2) Model 4 RR maps (Temp, Rain, Flood, Humidity) with improved readable legends
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
# FIXED: Improved error handling when searching for parameter matches
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

computeRR <- function(net_effect_matrix) {
  data.frame(RR = exp(colMeans(net_effect_matrix, na.rm = TRUE)))
}

# Enhanced color palette
red_palette <- colorRampPalette(
  c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
)(200)

blue_palette <- colorRampPalette(
  c("#F1EEF6", "#BDC9E1", "#74A9CF", "#2B8CBE", "#045A8D")
)(200)

base_map_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.background  = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA)
  )

# Format numbers to avoid scientific notation and handle different scales
format_numbers <- function(x) {
  if(any(x >= 1000000, na.rm=TRUE)) {
    # Use K for thousands and M for millions
    sapply(x, function(n) {
      if(n >= 1000000) {
        paste0(round(n/1000000, 1), "M")
      } else if(n >= 1000) {
        paste0(round(n/1000, 1), "K")
      } else {
        format(n, digits=2, scientific=FALSE)
      }
    })
  } else if(any(x >= 10, na.rm=TRUE)) {
    # For values above 10, round to 1 decimal place
    format(round(x, 1), scientific=FALSE)
  } else {
    # For small values, keep 2 decimal places
    format(round(x, 2), scientific=FALSE)
  }
}

# Map plotting function with enhanced legend control
plot_map <- function(
    sf_object, var_col, fill_label,
    title_str    = NULL,
    subtitle_str = NULL,
    fill_limits  = NULL,
    breaks       = NULL, 
    labels       = NULL,
    palette      = red_palette
) {
  # Format the labels to avoid scientific notation
  if(is.null(labels) && !is.null(breaks)) {
    labels <- format_numbers(breaks)
  }
  
  p <- ggplot(sf_object) +
    geom_sf(aes_string(fill = var_col), color="#666666", size=0.1) +
    scale_fill_gradientn(
      colors = palette,
      limits = fill_limits,
      oob    = scales::squish,
      name   = fill_label,
      breaks = breaks,
      labels = labels
    ) +
    labs(title=title_str, subtitle=subtitle_str) +
    base_map_theme +
    theme(
      plot.title        = element_text(size=14, face="bold", margin=margin(b=3)),
      plot.subtitle     = element_text(size=12, color="gray40", margin=margin(b=5)),
      legend.key.width  = unit(1.0, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.title      = element_text(size=10, face="bold"),
      legend.text       = element_text(size=9),
      legend.position   = "bottom",
      legend.margin     = margin(t=0, r=0, b=5, l=0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.x  = unit(0.2, "cm"),
      plot.margin       = margin(t=5, r=15, b=15, l=15)
    )
  
  p
}

#-----------------------------
# 4) Summaries
#-----------------------------
n_regions <- length(unique(shp_sf$LHD_code))

# Check for valid data before processing
check_valid_data <- function(net_effect, variable_name) {
  valid_pct <- sum(!is.na(net_effect)) / length(net_effect) * 100
  cat(sprintf("%s has %.1f%% valid data\n", variable_name, valid_pct))
  
  if(valid_pct < 1) {
    warning(sprintf("Very few valid values for %s (%.1f%%). Check model parameters.", 
                    variable_name, valid_pct))
  }
  return(valid_pct > 0)  # Return TRUE if we have at least some valid data
}

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

m3_df <- data.frame(
  LHD_code  = 1:n_regions,
  RR_temp   = if(valid_m3_temp) computeRR(m3_temp_net)$RR else rep(NA, n_regions),
  RR_rain   = if(valid_m3_rain) computeRR(m3_rain_net)$RR else rep(NA, n_regions),
  RR_flood  = if(valid_m3_flood) computeRR(m3_flood_net)$RR else rep(NA, n_regions),
  RR_humid  = if(valid_m3_humid) computeRR(m3_humid_net)$RR else rep(NA, n_regions)
)

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

m4_df <- data.frame(
  LHD_code  = 1:n_regions,
  RR_temp   = if(valid_m4_temp) computeRR(m4_temp_net)$RR else rep(NA, n_regions),
  RR_rain   = if(valid_m4_rain) computeRR(m4_rain_net)$RR else rep(NA, n_regions),
  RR_flood  = if(valid_m4_flood) computeRR(m4_flood_net)$RR else rep(NA, n_regions),
  RR_humid  = if(valid_m4_humid) computeRR(m4_humid_net)$RR else rep(NA, n_regions)
)

# Apply reasonable limits to handle outliers
cap_values <- function(x, cap_low = NULL, cap_high = NULL) {
  if(all(is.na(x))) return(x) # Return as is if all values are NA
  if(!is.null(cap_low)) x <- pmax(x, cap_low, na.rm = TRUE)
  if(!is.null(cap_high)) x <- pmin(x, cap_high, na.rm = TRUE)
  return(x)
}

# Check outliers and potentially cap extreme values
outlier_check <- function(x, title="Value") {
  if(all(is.na(x))) {
    cat(sprintf("%s - All values are NA\n", title))
    return(c(lower=NA, upper=NA))
  }
  
  cat(sprintf("%s - Range: [%.2f, %.2f], Mean: %.2f, Median: %.2f\n", 
              title, min(x, na.rm=TRUE), max(x, na.rm=TRUE), 
              mean(x, na.rm=TRUE), median(x, na.rm=TRUE)))
  
  # Check for potential outliers
  q1 <- quantile(x, 0.25, na.rm=TRUE)
  q3 <- quantile(x, 0.75, na.rm=TRUE)
  iqr <- q3 - q1
  upper_bound <- q3 + 1.5 * iqr
  lower_bound <- q1 - 1.5 * iqr
  
  outliers <- x[(x > upper_bound | x < lower_bound) & !is.na(x)]
  if(length(outliers) > 0) {
    cat(sprintf("  Found %d potential outliers outside [%.2f, %.2f]\n", 
                length(outliers), lower_bound, upper_bound))
  }
  
  return(c(lower=lower_bound, upper=upper_bound))
}

# Display some summary information
cat("\nModel 3 RR Summaries:\n")
m3_temp_bounds  <- outlier_check(m3_df$RR_temp, "Temperature RR")
m3_rain_bounds  <- outlier_check(m3_df$RR_rain, "Rainfall RR") 
m3_flood_bounds <- outlier_check(m3_df$RR_flood, "Flood RR")
m3_humid_bounds <- outlier_check(m3_df$RR_humid, "Humidity RR")

cat("\nModel 4 RR Summaries:\n")
m4_temp_bounds  <- outlier_check(m4_df$RR_temp, "Temperature RR")
m4_rain_bounds  <- outlier_check(m4_df$RR_rain, "Rainfall RR")
m4_flood_bounds <- outlier_check(m4_df$RR_flood, "Flood RR")
m4_humid_bounds <- outlier_check(m4_df$RR_humid, "Humidity RR")

# Optionally cap extreme values for better visualization
# Note: This is only for visualization - the actual values remain in the dataframe
m3_df$RR_temp_capped  <- cap_values(m3_df$RR_temp, m3_temp_bounds["lower"], m3_temp_bounds["upper"])
m3_df$RR_rain_capped  <- cap_values(m3_df$RR_rain, m3_rain_bounds["lower"], m3_rain_bounds["upper"])
m3_df$RR_flood_capped <- cap_values(m3_df$RR_flood, m3_flood_bounds["lower"], m3_flood_bounds["upper"])
m3_df$RR_humid_capped <- cap_values(m3_df$RR_humid, m3_humid_bounds["lower"], m3_humid_bounds["upper"])

m4_df$RR_temp_capped  <- cap_values(m4_df$RR_temp, m4_temp_bounds["lower"], m4_temp_bounds["upper"])
m4_df$RR_rain_capped  <- cap_values(m4_df$RR_rain, m4_rain_bounds["lower"], m4_rain_bounds["upper"])
m4_df$RR_flood_capped <- cap_values(m4_df$RR_flood, m4_flood_bounds["lower"], m4_flood_bounds["upper"])
m4_df$RR_humid_capped <- cap_values(m4_df$RR_humid, m4_humid_bounds["lower"], m4_humid_bounds["upper"])

# Join with shapefile
m3_map <- left_join(shp_sf, m3_df, by="LHD_code")
m4_map <- left_join(shp_sf, m4_df, by="LHD_code")

#-----------------------------
# 5) MODEL 3: Relative Risk Maps
#-----------------------------
# Calculate ranges for each variable
# Use either capped or uncapped values
use_capped <- TRUE  # Set to FALSE to use uncapped values

if(use_capped) {
  temp_range  <- range(m3_map$RR_temp_capped,  na.rm=TRUE)
  rain_range  <- range(m3_map$RR_rain_capped,  na.rm=TRUE)
  flood_range <- range(m3_map$RR_flood_capped, na.rm=TRUE)
  humid_range <- range(m3_map$RR_humid_capped, na.rm=TRUE)
} else {
  temp_range  <- range(m3_map$RR_temp,  na.rm=TRUE)
  rain_range  <- range(m3_map$RR_rain,  na.rm=TRUE)
  flood_range <- range(m3_map$RR_flood, na.rm=TRUE)
  humid_range <- range(m3_map$RR_humid, na.rm=TRUE)
}

# Define default ranges if data is entirely missing
set_default_range <- function(range_vals, default_min=0.5, default_max=2) {
  if(all(is.na(range_vals)) || range_vals[1] > range_vals[2]) {
    return(c(default_min, default_max))
  }
  return(range_vals)
}

# Define exactly 5 breaks for each variable with proper scaling
create_breaks <- function(range_vals, n_breaks=5, default_min=0.5, default_max=2) {
  # If range is invalid, use defaults
  if(all(is.na(range_vals)) || range_vals[1] > range_vals[2]) {
    range_vals <- c(default_min, default_max)
  }
  
  # Special handling for extreme values
  if(max(range_vals, na.rm=TRUE) > 1000) {
    # Round to nearest million or thousand for readability
    magnitude <- 10^floor(log10(max(range_vals, na.rm=TRUE)))
    return(seq(
      from = floor(range_vals[1]/magnitude) * magnitude,
      to = ceiling(range_vals[2]/magnitude) * magnitude,
      length.out = n_breaks
    ))
  } else {
    # Use regular sequence with sensible rounding
    return(seq(
      from = floor(range_vals[1] * 10) / 10,
      to = ceiling(range_vals[2] * 10) / 10,
      length.out = n_breaks
    ))
  }
}

temp_breaks  <- create_breaks(temp_range)
rain_breaks  <- create_breaks(rain_range)
flood_breaks <- create_breaks(flood_range)
humid_breaks <- create_breaks(humid_range)

# Create plots using either capped or uncapped values - with NA checking
m3_temp_plot <- if(valid_m3_temp) {
  plot_map(
    m3_map, 
    ifelse(use_capped, "RR_temp_capped", "RR_temp"), 
    "Relative Risk",
    title_str   = "Temperature",
    subtitle_str = "Model 3",
    fill_limits = temp_range,
    breaks      = temp_breaks
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid temperature data") +
    theme_void() +
    labs(title="Temperature", subtitle="Model 3")
}

m3_rain_plot <- if(valid_m3_rain) {
  plot_map(
    m3_map, 
    ifelse(use_capped, "RR_rain_capped", "RR_rain"), 
    "Relative Risk",
    title_str   = "Rainfall",
    subtitle_str = "Model 3",
    fill_limits = rain_range,
    breaks      = rain_breaks,
    palette     = red_palette
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid rainfall data") +
    theme_void() +
    labs(title="Rainfall", subtitle="Model 3")
}

m3_flood_plot <- if(valid_m3_flood) {
  plot_map(
    m3_map, 
    ifelse(use_capped, "RR_flood_capped", "RR_flood"), 
    "Relative Risk",
    title_str   = "Flooding",
    subtitle_str = "Model 3",
    fill_limits = flood_range,
    breaks      = flood_breaks
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid flood data") +
    theme_void() +
    labs(title="Flooding", subtitle="Model 3")
}

m3_humid_plot <- if(valid_m3_humid) {
  plot_map(
    m3_map, 
    ifelse(use_capped, "RR_humid_capped", "RR_humid"), 
    "Relative Risk",
    title_str   = "Humidity",
    subtitle_str = "Model 3",
    fill_limits = humid_range,
    breaks      = humid_breaks,
    palette     = red_palette
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid humidity data") +
    theme_void() +
    labs(title="Humidity", subtitle="Model 3")
}

# Create first page with 2×2 layout - handling NA cases
if(sum(c(valid_m3_temp, valid_m3_rain, valid_m3_flood, valid_m3_humid)) > 0) {
  # Create a list to hold valid plots
  m3_plots <- list()
  
  # Only add plots for variables with valid data
  if(valid_m3_temp) m3_plots$temp <- m3_temp_plot
  if(valid_m3_rain) m3_plots$rain <- m3_rain_plot
  if(valid_m3_flood) m3_plots$flood <- m3_flood_plot
  if(valid_m3_humid) m3_plots$humid <- m3_humid_plot
  
  # Set up the layout depending on how many plots we have
  if(length(m3_plots) == 4) {
    m3_rr_page1 <- (m3_plots$temp + m3_plots$rain) / (m3_plots$flood + m3_plots$humid)
  } else if(length(m3_plots) == 3) {
    # Create a 3-plot layout
    plot_names <- names(m3_plots)
    m3_rr_page1 <- m3_plots[[plot_names[1]]] + m3_plots[[plot_names[2]]] + m3_plots[[plot_names[3]]]
  } else if(length(m3_plots) == 2) {
    m3_rr_page1 <- m3_plots[[1]] + m3_plots[[2]]
  } else if(length(m3_plots) == 1) {
    m3_rr_page1 <- m3_plots[[1]]
  } else {
    # Create a placeholder plot
    m3_rr_page1 <- ggplot() + 
      geom_text(aes(x=0, y=0, label="No valid data for Model 3")) +
      theme_void()
  }
  
  # Add annotation
  m3_rr_page1 <- m3_rr_page1 + 
    plot_annotation(
      title="Model 3: Salmonella Relative Risk by Exposure Type",
      theme=theme(
        plot.title=element_text(size=18, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=14, hjust=0.5),
        plot.margin=margin(t=10, b=10)
      )
    ) &
    theme(
      plot.margin = margin(5, 5, 20, 5)
    )
} else {
  # Create a placeholder if no valid data
  m3_rr_page1 <- ggplot() + 
    geom_text(aes(x=0, y=0, label="No valid data for Model 3")) +
    theme_void() +
    theme(
      plot.title=element_text(size=18, face="bold", hjust=0.5),
      plot.margin=margin(t=10, b=10)
    )
}

#-----------------------------
# 6) MODEL 4: Relative Risk Maps
#-----------------------------
# Calculate ranges for each variable
if(use_capped) {
  temp_range4  <- set_default_range(range(m4_map$RR_temp_capped,  na.rm=TRUE))
  rain_range4  <- set_default_range(range(m4_map$RR_rain_capped,  na.rm=TRUE))
  flood_range4 <- set_default_range(range(m4_map$RR_flood_capped, na.rm=TRUE))
  humid_range4 <- set_default_range(range(m4_map$RR_humid_capped, na.rm=TRUE))
} else {
  temp_range4  <- set_default_range(range(m4_map$RR_temp,  na.rm=TRUE))
  rain_range4  <- set_default_range(range(m4_map$RR_rain,  na.rm=TRUE))
  flood_range4 <- set_default_range(range(m4_map$RR_flood, na.rm=TRUE))
  humid_range4 <- set_default_range(range(m4_map$RR_humid, na.rm=TRUE))
}

# Define exactly 5 breaks for each variable with proper scaling
temp_breaks4  <- create_breaks(temp_range4)
rain_breaks4  <- create_breaks(rain_range4)
flood_breaks4 <- create_breaks(flood_range4)
humid_breaks4 <- create_breaks(humid_range4)

# Create plots for Model 4 - with NA checking
m4_temp_plot <- if(valid_m4_temp) {
  plot_map(
    m4_map, 
    ifelse(use_capped, "RR_temp_capped", "RR_temp"), 
    "Relative Risk",
    title_str   = "Temperature",
    subtitle_str = "Model 4",
    fill_limits = temp_range4,
    breaks      = temp_breaks4
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid temperature data") +
    theme_void() +
    labs(title="Temperature", subtitle="Model 4")
}

m4_rain_plot <- if(valid_m4_rain) {
  plot_map(
    m4_map, 
    ifelse(use_capped, "RR_rain_capped", "RR_rain"), 
    "Relative Risk",
    title_str   = "Rainfall",
    subtitle_str = "Model 4",
    fill_limits = rain_range4,
    breaks      = rain_breaks4,
    palette     = red_palette
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid rainfall data") +
    theme_void() +
    labs(title="Rainfall", subtitle="Model 4")
}

m4_flood_plot <- if(valid_m4_flood) {
  plot_map(
    m4_map, 
    ifelse(use_capped, "RR_flood_capped", "RR_flood"), 
    "Relative Risk",
    title_str   = "Flooding",
    subtitle_str = "Model 4",
    fill_limits = flood_range4,
    breaks      = flood_breaks4
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid flood data") +
    theme_void() +
    labs(title="Flooding", subtitle="Model 4")
}

m4_humid_plot <- if(valid_m4_humid) {
  plot_map(
    m4_map, 
    ifelse(use_capped, "RR_humid_capped", "RR_humid"), 
    "Relative Risk",
    title_str   = "Humidity",
    subtitle_str = "Model 4",
    fill_limits = humid_range4,
    breaks      = humid_breaks4,
    palette     = red_palette
  )
} else {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label="No valid humidity data") +
    theme_void() +
    labs(title="Humidity", subtitle="Model 4")
}

# Create second page with 2×2 layout - handling NA cases
if(sum(c(valid_m4_temp, valid_m4_rain, valid_m4_flood, valid_m4_humid)) > 0) {
  # Create a list to hold valid plots
  m4_plots <- list()
  
  # Only add plots for variables with valid data
  if(valid_m4_temp) m4_plots$temp <- m4_temp_plot
  if(valid_m4_rain) m4_plots$rain <- m4_rain_plot
  if(valid_m4_flood) m4_plots$flood <- m4_flood_plot
  if(valid_m4_humid) m4_plots$humid <- m4_humid_plot
  
  # Set up the layout depending on how many plots we have
  if(length(m4_plots) == 4) {
    m4_rr_page1 <- (m4_plots$temp + m4_plots$rain) / (m4_plots$flood + m4_plots$humid)
  } else if(length(m4_plots) == 3) {
    # Create a 3-plot layout
    plot_names <- names(m4_plots)
    m4_rr_page1 <- m4_plots[[plot_names[1]]] + m4_plots[[plot_names[2]]] + m4_plots[[plot_names[3]]]
  } else if(length(m4_plots) == 2) {
    m4_rr_page1 <- m4_plots[[1]] + m4_plots[[2]]
  } else if(length(m4_plots) == 1) {
    m4_rr_page1 <- m4_plots[[1]]
  } else {
    # Create a placeholder plot
    m4_rr_page1 <- ggplot() + 
      geom_text(aes(x=0, y=0, label="No valid data for Model 4")) +
      theme_void()
  }
  
  # Add annotation
  m4_rr_page1 <- m4_rr_page1 + 
    plot_annotation(
      title="Model 4: Salmonella Relative Risk by Exposure Type",
      theme=theme(
        plot.title=element_text(size=18, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=14, hjust=0.5),
        plot.margin=margin(t=10, b=10)
      )
    ) &
    theme(
      plot.margin = margin(5, 5, 20, 5)
    )
} else {
  # Create a placeholder if no valid data
  m4_rr_page1 <- ggplot() + 
    geom_text(aes(x=0, y=0, label="No valid data for Model 4")) +
    theme_void() +
    theme(
      plot.title=element_text(size=18, face="bold", hjust=0.5),
      plot.margin=margin(t=10, b=10)
    )
}

#-----------------------------
# 7) Output to PDF - just the RR maps
#-----------------------------
pdf("model_output/relative_risk_maps_nimble.pdf", width=14, height=10)

# Page 1: Model 3 RR (2×2)
print(m3_rr_page1)

# Page 2: Model 4 RR (2×2)
print(m4_rr_page1)

dev.off()

cat("\nRelative Risk maps for NIMBLE models saved to 'model_output/relative_risk_maps_nimble.pdf'.\n")

# Also save as separate images for potential inclusion in documents
dir.create("model_output/risk_maps", showWarnings = FALSE)

# Save individual plots as PNG, only if they contain valid data
if(valid_m3_temp) ggsave("model_output/risk_maps/model3_temperature.png", m3_temp_plot, width=7, height=6, dpi=300)
if(valid_m3_rain) ggsave("model_output/risk_maps/model3_rainfall.png", m3_rain_plot, width=7, height=6, dpi=300)
if(valid_m3_flood) ggsave("model_output/risk_maps/model3_flooding.png", m3_flood_plot, width=7, height=6, dpi=300)
if(valid_m3_humid) ggsave("model_output/risk_maps/model3_humidity.png", m3_humid_plot, width=7, height=6, dpi=300)

if(valid_m4_temp) ggsave("model_output/risk_maps/model4_temperature.png", m4_temp_plot, width=7, height=6, dpi=300)
if(valid_m4_rain) ggsave("model_output/risk_maps/model4_rainfall.png", m4_rain_plot, width=7, height=6, dpi=300)
if(valid_m4_flood) ggsave("model_output/risk_maps/model4_flooding.png", m4_flood_plot, width=7, height=6, dpi=300)
if(valid_m4_humid) ggsave("model_output/risk_maps/model4_humidity.png", m4_humid_plot, width=7, height=6, dpi=300)

cat("Individual map images saved to 'model_output/risk_maps/' folder.\n")