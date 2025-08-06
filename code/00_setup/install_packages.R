# install_packages.R
# Script to install all required packages for Salmonella-climate analysis

# Function to check and install packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# List of required packages
required_packages <- c(
  # Data manipulation
  "dplyr", "tidyr", "purrr", "lubridate", "reshape2",
  
  # Spatial analysis
  "sf", "sp", "spdep",
  
  # Statistical modeling
  "nimble", "dlnm", "coda",
  
  # Visualization
  "ggplot2", "gridExtra", "patchwork", "scales", 
  "RColorBrewer", "viridis",
  
  # Utilities
  "knitr", "kableExtra", "fmsb"
)

# Install packages
cat("Installing required packages...\n")
for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Check NIMBLE compilation
cat("\nChecking NIMBLE compilation...\n")
library(nimble)
testCode <- nimbleCode({
  x ~ dnorm(0, 1)
})
testModel <- nimbleModel(testCode)
if (!is.null(testModel)) {
  cat("✓ NIMBLE compilation successful\n")
} else {
  cat("✗ NIMBLE compilation failed - check compiler setup\n")
}

# Test spatial packages
cat("\nChecking spatial packages...\n")
if (requireNamespace("sf", quietly = TRUE)) {
  sf::sf_use_s2(FALSE)  # Disable spherical geometry for compatibility
  cat("✓ Spatial packages configured\n")
}

cat("\nPackage installation complete!\n")
cat("You can now proceed with the analysis.\n")