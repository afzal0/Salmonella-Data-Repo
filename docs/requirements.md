# Software Requirements and Dependencies

## R Environment

### R Version
- **Minimum**: R 4.2.0
- **Recommended**: R 4.3.0 or higher
- **Tested on**: R 4.3.2

### RStudio
- **Recommended**: RStudio 2023.06.0 or higher
- **Optional but helpful** for interactive analysis

## Required R Packages

### Core Data Manipulation
```R
install.packages(c(
  "dplyr",      # >= 1.1.0
  "tidyr",      # >= 1.3.0
  "purrr",      # >= 1.0.0
  "lubridate",  # >= 1.9.0
  "reshape2"    # >= 1.4.4
))
```

### Spatial Analysis
```R
install.packages(c(
  "sf",         # >= 1.0-12
  "sp",         # >= 1.6-0
  "spdep"       # >= 1.2-8
))
```

### Statistical Modeling
```R
# NIMBLE requires Rtools on Windows
install.packages(c(
  "nimble",     # >= 1.0.0
  "dlnm",       # >= 2.4.7
  "splines",    # (base R)
  "coda"        # >= 0.19-4
))
```

### Visualization
```R
install.packages(c(
  "ggplot2",      # >= 3.4.0
  "gridExtra",    # >= 2.3
  "patchwork",    # >= 1.1.2
  "scales",       # >= 1.2.1
  "RColorBrewer", # >= 1.1-3
  "viridis"       # >= 0.6.2
))
```

### Additional Utilities
```R
install.packages(c(
  "knitr",      # >= 1.42
  "kableExtra", # >= 1.3.4
  "fmsb"        # >= 0.7.5
))
```

## System Requirements

### Operating System
- **Supported**: Windows 10+, macOS 10.15+, Ubuntu 20.04+
- **Tested on**: Windows 11, macOS Ventura, Ubuntu 22.04

### Hardware Requirements
- **CPU**: Multi-core processor recommended (for MCMC parallelization)
- **RAM**: 
  - Minimum: 8 GB
  - Recommended: 16 GB or more
- **Storage**: At least 5 GB free space for:
  - Data files (~500 MB)
  - Model outputs (~1.5 GB)
  - Temporary files during computation

### Compiler Requirements (for NIMBLE)
- **Windows**: Rtools43 or appropriate version for your R
- **macOS**: Xcode Command Line Tools
- **Linux**: gcc/g++ compiler suite

## Installation Script

Save and run this script to install all dependencies:

```R
# install_packages.R

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

cat("\nPackage installation complete!\n")
```

## Environment Setup

### Setting Memory Limits (if needed)
```R
# Increase memory limit on Windows
memory.limit(size = 16000)  # 16GB

# For Unix-based systems, set before starting R:
# export R_MAX_VSIZE=16g
```

### Parallel Processing Setup
```R
# Detect available cores
n_cores <- parallel::detectCores() - 1

# Configure NIMBLE for parallel MCMC
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
nimbleOptions(MCMCorderPosteriorPredictiveSamplers = FALSE)
```

## Troubleshooting Installation

### NIMBLE Installation Issues

#### Windows
1. Install Rtools from: https://cran.r-project.org/bin/windows/Rtools/
2. Add Rtools to PATH
3. Restart R and retry installation

#### macOS
1. Install Xcode Command Line Tools:
   ```bash
   xcode-select --install
   ```
2. If issues persist, install full Xcode from App Store

#### Linux
1. Install build essentials:
   ```bash
   sudo apt-get update
   sudo apt-get install build-essential
   ```

### Spatial Package Issues

If `sf` installation fails:
```bash
# Ubuntu/Debian
sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev

# macOS (with Homebrew)
brew install gdal proj geos

# Then retry in R
install.packages("sf")
```

## Version Control

### Package Version Locking
To ensure reproducibility, save package versions:

```R
# Save current package versions
pkg_versions <- installed.packages()[required_packages, "Version"]
saveRDS(pkg_versions, "package_versions.rds")

# To restore specific versions later:
# devtools::install_version("package", version = "x.y.z")
```

### Using renv (Recommended)
```R
# Initialize renv for the project
install.packages("renv")
renv::init()

# Snapshot current environment
renv::snapshot()

# Others can restore with:
# renv::restore()
```

## Docker Alternative

For complete reproducibility, a Docker image is available:

```dockerfile
# Dockerfile
FROM rocker/geospatial:4.3.2

RUN install2.r --error \
    nimble dlnm coda patchwork fmsb kableExtra

WORKDIR /analysis
```

Build and run:
```bash
docker build -t salmonella-analysis .
docker run -it -v $(pwd):/analysis salmonella-analysis
```

## Contact for Technical Support

For installation or technical issues:
- GitHub Issues: [repository]/issues
- Email: [technical contact]