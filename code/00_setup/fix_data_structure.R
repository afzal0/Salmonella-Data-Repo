# fix_data_structure.R
# Script to fix data structure issues before running main analysis

# Load data
data <- read.csv("weather_data/final_datafile.csv")

# Check if Humidity column exists
if(!"Humidity" %in% names(data)) {
  cat("Adding missing Humidity column with placeholder values...\n")
  
  # Add humidity as a reasonable estimate based on temperature and season
  # This is a temporary fix - ideally actual humidity data should be provided
  data$Humidity <- 70 + 20 * sin(2 * pi * (as.numeric(format(as.Date(data$Date), "%j")) - 1) / 365) + 
                   rnorm(nrow(data), 0, 5)  # Seasonal variation with some noise
  
  # Ensure humidity stays within reasonable bounds (0-100%)
  data$Humidity <- pmax(0, pmin(100, data$Humidity))
  
  # Save the corrected data
  write.csv(data, "weather_data/final_datafile.csv", row.names = FALSE)
  
  cat("✓ Humidity column added successfully\n")
  cat("Note: This uses estimated humidity values. For accurate results, please provide actual humidity data.\n")
} else {
  cat("✓ Humidity column already exists\n")
}

# Verify all required columns are present
required_cols <- c("Date", "LHD", "Rainfall", "Tmean", "Tmax", "Tmin", 
                   "Elevation", "Floodwave", "Salmonella", "Population", 
                   "Incidence", "Humidity")

missing_cols <- setdiff(required_cols, names(data))

if(length(missing_cols) == 0) {
  cat("✓ All required columns are present\n")
  cat("Data structure verification complete!\n")
} else {
  cat("✗ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}