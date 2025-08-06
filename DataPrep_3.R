###############################################################
#   DATA PREPARATION SCRIPT FOR SPATIAL B-DLNM (Models 3 & 4)
###############################################################

# 1) Load Required Libraries
library(dplyr)
library(tidyr)
library(lubridate)   # For parsing dates
library(sf)          # For reading/writing spatial data
library(sp)          # For Spatial* objects
library(spdep)       # For creating neighbor/adjacency structures
library(dlnm)        # For distributed lag non-linear model crossbasis
library(splines)     # For spline functions

# Ensure "model_input" directory exists
dir.create("model_input", showWarnings = FALSE)

#---------------------------------------------------------------------
# 2) Read the Epidemiological CSV Data
#---------------------------------------------------------------------
# This CSV file is expected to have columns like:
# Date, LHD, Rainfall, Tmean, Tmax, Tmin, Elevation, Floodwave,
# Salmonella, Population, Incidence, Humidity
my_data <- read.csv("weather_data/final_datafile.csv", stringsAsFactors = FALSE)

# Convert Date to a proper Date class
my_data$Date <- as.Date(my_data$Date, format = "%Y-%m-%d")

# Inspect basic structure (optional)
str(my_data)

# Optional: check for missing data
cat("Missing values per column:\n")
print(colSums(is.na(my_data)))

#---------------------------------------------------------------------
# 3) Sort & Create a Numeric LHD Code
#---------------------------------------------------------------------
my_data <- my_data %>%
  arrange(LHD, Date) %>%
  mutate(
    LHD_code = as.integer(factor(LHD))  # factor -> integer
  )

cat("\nUnique LHDs in data:\n")
print(unique(my_data$LHD))

#---------------------------------------------------------------------
# 4) Create TWO Data Frames:
#    a) data_cco for Case-Crossover (Model 3)
#    b) data_ts  for Time-Series   (Model 4)
#---------------------------------------------------------------------

### (a) data_cco
data_cco <- my_data %>%
  mutate(
    strata = paste(
      LHD_code,
      format(Date, "%Y-%m"),
      sep = ":"
    )
  )

### (b) data_ts
data_ts <- my_data %>%
  mutate(
    month = month(Date),
    year  = year(Date)
  )

#---------------------------------------------------------------------
# 5) Load Shapefile and Build Adjacency (Spatial Component)
#---------------------------------------------------------------------
# Example shapefile: "NSW_LHD_Boundaries.shp"
shp <- st_read("salmonella_hassan/spatial_data/NSW_LHD_Boundaries.shp")

# Create an LHD column in the shapefile to match the data
shp$LHD <- as.character(shp$lhd_name)

# Convert to sp
shp_sp <- as(shp, "Spatial")

# Make a numeric code in shapefile that aligns with data's LHD
shp_sp@data$LHD_code <- as.integer(factor(shp_sp@data$LHD))

# Create neighbor/adjacency structures
nb_list    <- poly2nb(shp_sp, row.names = shp_sp@data$LHD_code)
nb_weights <- nb2listw(nb_list, style = "B")  # binary weights
adj_mat    <- nb2mat(nb_list, style = "B", zero.policy = TRUE)

# Store adjacency in a form suitable for iCAR in JAGS:
nb_info <- nb2WB(nb_list)  # yields: $adj, $weights, $num
adj_vec <- nb_info$adj
wts_vec <- nb_info$weights
num_vec <- nb_info$num

# sumNumNeigh often used to constrain the iCAR (sum(num_vec)):
sumNumNeigh <- sum(num_vec)

#---------------------------------------------------------------------
# 6) CREATE CROSS-BASIS FOR DLNM VARIABLES
#    Tmean, Rainfall, Floodwave, and Humidity for both data_cco & data_ts
#---------------------------------------------------------------------
# Set a lag of 2 for all parameters
max_lag <- 2
var_knots <- c(0.10, 0.50, 0.90)

### Cross-Basis for data_cco

## Tmean
cb_cco_Tmean <- crossbasis(
  x   = data_cco$Tmean,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_cco$Tmean, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Rainfall
cb_cco_Rain <- crossbasis(
  x   = data_cco$Rainfall,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_cco$Rainfall, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Floodwave
cb_cco_Flood <- crossbasis(
  x   = data_cco$Floodwave,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_cco$Floodwave, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Humidity
cb_cco_Humidity <- crossbasis(
  x   = data_cco$Humidity,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_cco$Humidity, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

# Append the cross-basis expansions to data_cco
data_cco <- cbind(
  data_cco,
  as.data.frame(cb_cco_Tmean),
  as.data.frame(cb_cco_Rain),
  as.data.frame(cb_cco_Flood),
  as.data.frame(cb_cco_Humidity)
)

### Cross-Basis for data_ts

## Tmean
cb_ts_Tmean <- crossbasis(
  x   = data_ts$Tmean,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_ts$Tmean, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Rainfall
cb_ts_Rain <- crossbasis(
  x   = data_ts$Rainfall,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_ts$Rainfall, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Floodwave
cb_ts_Flood <- crossbasis(
  x   = data_ts$Floodwave,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_ts$Floodwave, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

## Humidity
cb_ts_Humidity <- crossbasis(
  x   = data_ts$Humidity,
  lag = max_lag,
  argvar = list(fun = "ns", knots = quantile(data_ts$Humidity, var_knots, na.rm = TRUE)),
  arglag = list(fun = "ns")
)

# Append the cross-basis expansions to data_ts
data_ts <- cbind(
  data_ts,
  as.data.frame(cb_ts_Tmean),
  as.data.frame(cb_ts_Rain),
  as.data.frame(cb_ts_Flood),
  as.data.frame(cb_ts_Humidity)
)

#---------------------------------------------------------------------
# 7) SAVE ALL PREPARED OBJECTS FOR MODEL 3 & 4
#    in the "model_input" directory
#---------------------------------------------------------------------

# (a) Data frames
save(data_cco, file = "model_input/data_casecrossover.RData")  # For Model 3
save(data_ts,  file = "model_input/data_timeseries.RData")     # For Model 4

# (b) Spatial adjacency objects
save(
  shp_sp, nb_list, nb_weights, adj_mat,
  nb_info, adj_vec, wts_vec, num_vec, sumNumNeigh,
  file = "model_input/spatial_adjacency.RData"
)

# Optionally, also save the cross-basis objects separately if desired:
# save(cb_cco_Tmean, cb_cco_Rain, cb_cco_Flood, cb_cco_Humidity,
#      cb_ts_Tmean, cb_ts_Rain, cb_ts_Flood, cb_ts_Humidity,
#      file = "model_input/crossbasis_objects.RData")

cat("\n------------------------------------------------\n")
cat("DATA PREPARATION COMPLETE (Models 3 & 4).\n\n")
cat("Saved in 'model_input' folder:\n")
cat(" 1) data_casecrossover.RData  (data_cco)\n")
cat(" 2) data_timeseries.RData     (data_ts)\n")
cat(" 3) spatial_adjacency.RData   (spatial adjacency objects)\n")
cat("\nUse these files in your SB-DLNM scripts (Model 3 & 4).\n")
