# Code Documentation

This directory contains all R scripts for replicating the analysis of climate impacts on Salmonella infections in NSW.

## Directory Structure

- `00_setup/`: Package installation and environment setup
- `01_data_preparation/`: Data cleaning and preprocessing
- `02_modeling/`: Statistical model implementation
- `03_visualization/`: Results visualization and mapping

## Analysis Pipeline

The scripts should be executed in the following order:

### 1. Setup Environment
```R
source("code/00_setup/install_packages.R")
```

### 2. Data Preparation
```R
source("code/01_data_preparation/DataPrep_3.R")
# Optional: Check data quality
source("code/01_data_preparation/Check_data.R")
```

### 3. Model Fitting
```R
source("code/02_modeling/SB_DLNM-New.R")
```
- Fits both Model 3 (case-crossover) and Model 4 (time-series)
- Runtime: ~2-3 hours depending on hardware
- Outputs saved to `results/models/`

### 4. Results Generation
```R
# Generate exposure-response curves
source("code/03_visualization/DLNM-Curves-2.R")

# Create spatial risk maps
source("code/03_visualization/Relative Risk Final.R")

# Generate probability exceedance maps
source("code/03_visualization/Probabibility Maps.R")
```

### 5. Model Evaluation
```R
source("code/02_modeling/Accuracy Assessment.R")
```

## Script Descriptions

### Data Preparation

#### DataPrep_3.R
- **Purpose**: Main data preparation script
- **Inputs**: 
  - `data/raw/datafinal.csv`
  - Shapefiles in `data/raw/`
- **Outputs**:
  - `data/processed/data_casecrossover.RData`
  - `data/processed/data_timeseries.RData`
  - `data/processed/spatial_adjacency.RData`
- **Key Steps**:
  1. Load and merge spatial/temporal data
  2. Create DLNM cross-basis matrices
  3. Generate spatial adjacency structures
  4. Format data for both models

#### Check_data.R
- **Purpose**: Data quality assessment
- **Inputs**: Processed RData files
- **Outputs**: Console diagnostics
- **Checks**: Missing values, duplicates, spatial structure integrity

### Modeling

#### SB_DLNM-New.R
- **Purpose**: Fit Bayesian hierarchical models
- **Inputs**: Processed data files
- **Outputs**:
  - `results/models/model3_results_nimble.RData`
  - `results/models/model4_results_nimble.RData`
- **Models**:
  - Model 3: Case-crossover with spatial effects
  - Model 4: Time-series with spatial effects
- **MCMC Settings**:
  - Chains: 3
  - Iterations: 5000
  - Burn-in: 1000
  - Thinning: 2

#### Accuracy Assessment.R
- **Purpose**: Evaluate model performance
- **Inputs**: Model results
- **Outputs**: 
  - Accuracy metrics (RMSE, R², MAPE, etc.)
  - Diagnostic plots
- **Assessments**:
  - Overall accuracy
  - Regional performance
  - Seasonal patterns
  - Case-load stratified metrics

### Visualization

#### DLNM-Curves-2.R
- **Purpose**: Generate exposure-response curves
- **Inputs**: Model results
- **Outputs**: `results/figures/climate_salmonella_associations.pdf`
- **Features**:
  - Separate curves by region and variable
  - 95% credible intervals
  - Reference lines at null effect

#### Relative Risk Final.R
- **Purpose**: Create spatial risk maps
- **Inputs**: Model results, shapefiles
- **Outputs**: 
  - `results/figures/relative_risk_maps_nimble.pdf`
  - Individual PNG files for each exposure
- **Map Types**:
  - Choropleth maps by exposure level
  - Consistent color scales across maps

#### Probabibility Maps.R
- **Purpose**: Map exceedance probabilities
- **Inputs**: Model results
- **Outputs**: `results/figures/relative_risk_probability_maps.pdf`
- **Thresholds**: P(RR > 1.0), P(RR > 1.2), P(RR > 1.5)

## Key Functions and Methods

### DLNM Cross-basis
```R
# Example cross-basis creation
cb_temp <- crossbasis(
  data$Tmean_centered,
  lag = 2,
  argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "ns", df = 3)
)
```

### NIMBLE Model Specification
```R
# Simplified model structure
model_code <- nimbleCode({
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha[strata[i]] + inprod(X[i,], beta[,region[i]]) + phi[region[i]]
  }
  
  # Spatial random effects
  phi[1:R] ~ dcar_normal(adj[], num[], tau.phi)
  
  # Priors
  # ...
})
```

### Relative Risk Calculation
```R
# Extract posterior samples
RR <- exp(beta_samples)
RR_mean <- apply(RR, 2:3, mean)
RR_CI <- apply(RR, 2:3, quantile, c(0.025, 0.975))
```

## Computational Requirements

- **RAM**: Minimum 8GB, recommended 16GB
- **Processing**: Multi-core recommended for MCMC
- **Storage**: ~2GB for all outputs
- **Time**: Full pipeline ~4-6 hours

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce MCMC iterations or thin more aggressively
2. **Convergence warnings**: Increase burn-in or iterations
3. **Missing packages**: Run install script again
4. **Spatial mismatches**: Check CRS alignment in shapefiles

### Debug Mode

Set `debug_mode <- TRUE` in scripts for:
- Verbose output
- Intermediate result saving
- Diagnostic plots

## Customization

### Modifying Models

To adjust model specifications:

1. Edit cross-basis parameters in `DataPrep_3.R`
2. Modify NIMBLE code in `SB_DLNM-New.R`
3. Update prior distributions as needed

### Adding New Visualizations

1. Create new script in `03_visualization/`
2. Follow naming convention
3. Document inputs/outputs
4. Update this README

## Contact

For code-related questions:
- Open an issue on GitHub
- Contact: [developer email]