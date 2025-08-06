# Detailed Methodology

## Study Design

### Overview
This study employs a spatial-temporal ecological design to investigate associations between climate variables and Salmonella incidence across 15 Local Health Districts (LHDs) in New South Wales, Australia, from January 1990 to December 2022.

### Study Area
- **Geographic scope**: New South Wales, Australia
- **Administrative units**: 15 Local Health Districts
- **Population**: Approximately 8.2 million (2022)
- **Area**: 800,642 km²

## Data Collection

### Epidemiological Data
- **Source**: NSW Notifiable Conditions Information Management System (NCIMS)
- **Case definition**: Laboratory-confirmed Salmonella infections
- **Temporal resolution**: Monthly aggregation
- **Spatial resolution**: Local Health District level
- **Period**: January 1991 - December 2022 (384 months)

### Climate Data
1. **Temperature**
   - Variables: Mean, maximum, and minimum monthly temperature
   - Source: Australian Bureau of Meteorology (BOM)
   - Processing: Station data interpolated to LHD centroids using inverse distance weighting

2. **Rainfall**
   - Variable: Total monthly rainfall
   - Source: BOM rain gauge network
   - Processing: Area-weighted average for each LHD

3. **Flood Events**
   - Variable: Flood severity index (0-10 scale)
   - Source: Emergency Management Australia disasters database
   - Classification: Based on affected area, duration, and impact

4. **Humidity**
   - Variable: Monthly average relative humidity
   - Source: BOM automatic weather stations
   - Processing: Spatial interpolation to LHD level

### Population Data
- **Source**: Australian Bureau of Statistics
- **Type**: Annual population estimates
- **Use**: Calculation of incidence rates (cases per 100,000)

## Statistical Analysis

### Model Framework

We implemented two complementary Bayesian hierarchical models using distributed lag non-linear models (DLNM) to capture complex exposure-response relationships and account for delayed effects.

### Model 3: Case-Crossover Design with Spatial Effects

#### Model Specification
```
Y_it | λ_it ~ Poisson(λ_it)
log(λ_it) = α_st + Σ_j Σ_l β_jl × cb_jl(X_jit) + φ_i
```

Where:
- Y_it: Salmonella count in LHD i at time t
- α_st: Stratum-specific intercept (LHD × Year-Month)
- cb_jl: Cross-basis function for exposure j at lag l
- φ_i: Spatial random effect for LHD i
- β_jl: Coefficient for exposure j at lag l

#### Key Features
- Controls for time-invariant confounders through stratification
- Captures within-stratum variation in exposures
- Spatial correlation via intrinsic CAR model

### Model 4: Time-Series with Spatial Effects

#### Model Specification
```
Y_it | λ_it ~ Poisson(λ_it)
log(λ_it) = α_im + Σ_j Σ_l β_jl × cb_jl(X_jit) + γ × trend_t + Σ_k δ_k × season_kt + φ_i
```

Additional terms:
- α_im: Region-month specific intercept
- γ: Linear trend coefficient
- δ_k: Seasonal harmonic coefficients
- trend_t: Linear time trend
- season_kt: Fourier seasonal terms

#### Key Features
- Explicit modeling of temporal trends
- Seasonal adjustment using Fourier terms
- Region-specific baseline risks by month

### Distributed Lag Non-Linear Models (DLNM)

#### Cross-basis Construction

For each climatic exposure, cross-basis functions were defined using natural cubic splines:
```R
crossbasis(
  x = exposure_centered,
  lag = 2,  # Maximum lag of 2 days
  argvar = list(fun = "ns", knots = quantile(x, c(0.1, 0.5, 0.9))),  # Knots at 10th, 50th, 90th percentiles
  arglag = list(fun = "ns", df = 3)   # Natural spline for lag dimension
)
```

#### Rationale
- Captures non-linear exposure-response relationships
- Models delayed effects up to 2 days (though aggregated to monthly level)
- Flexible functional forms via natural splines with knots at exposure percentiles
- All exposures mean-centered prior to basis construction

### Spatial Correlation

#### Intrinsic CAR Model
```
φ ~ ICAR(W, τ_φ)
```

Where:
- W: Spatial adjacency matrix (queen contiguity)
- τ_φ: Precision parameter
- ICAR: Intrinsic Conditional Autoregressive distribution

#### Adjacency Definition
- Queen contiguity: LHDs sharing any boundary point
- Binary weights: w_ij = 1 if adjacent, 0 otherwise
- Row-standardized for computation

### Prior Distributions

For both models, hierarchical normal distributions were used for regression coefficients:

```
# Hierarchical structure for exposure coefficients
β_r,j^(X) = μ_j^(X) + σ_j^(X) × θ_r,j^(X)
θ_r,j^(X) ~ N(0, τ_X^(-1))

# Hyperpriors
μ_j^(X) ~ N(0, 10³)
σ_j^(X) ~ U(0, 10)
τ_X ~ Gamma(0.1, 0.1)

# Other parameters
α ~ N(0, 10³)        # Intercepts
τ_φ ~ Gamma(0.1, 0.1)  # iCAR precision
```

### Model Implementation

#### MCMC Specifications
- **Software**: NIMBLE 1.0.0 in R
- **Chains**: 3 independent chains
- **Iterations**: 5,000 per chain
- **Burn-in**: 1,000 iterations
- **Thinning**: 2 (retain every 2nd sample)
- **Total posterior samples**: 6,000

#### Convergence Diagnostics
1. **Gelman-Rubin statistic**: R̂ < 1.1 for all parameters
2. **Effective sample size**: >1,000 for key parameters
3. **Trace plots**: Visual inspection for mixing
4. **Autocorrelation**: Check for excessive correlation

All models were fitted using NIMBLE package (version 1.0.0) in R version 4.3.3.

### Model Comparison

#### Deviance Information Criterion (DIC)
```
DIC = D̄ + pD
```
Where:
- D̄: Posterior mean deviance
- pD: Effective number of parameters

#### Quasi-AIC for Overdispersion
```
QAIC = -2L/ĉ + 2k
```
Where:
- L: Log-likelihood
- ĉ: Overdispersion parameter
- k: Number of parameters

### Risk Quantification

#### Relative Risk Calculation
```R
# For specific exposure level x at lag l
RR(x,l) = exp(Σ_j β_j × cb_j(x,l))

# Overall cumulative RR
RR_cumulative(x) = exp(Σ_l Σ_j β_jl × cb_jl(x))
```

#### Uncertainty Quantification
- 95% Credible Intervals from posterior distribution
- Exceedance probabilities: P(RR > threshold)
- Spatial uncertainty via posterior of φ

### Sensitivity Analyses

1. **Lag structure**: Test lags 0-3 months
2. **Degrees of freedom**: Vary df = 2-4 for splines
3. **Prior sensitivity**: Test alternative prior specifications
4. **Temporal stability**: Subset analysis by decade

## Software Implementation

### Computational Environment
- **Platform**: High-performance computing cluster
- **Parallelization**: 3 chains run in parallel
- **Memory**: 32GB RAM allocated per job
- **Runtime**: ~2-3 hours per model

### Reproducibility Features
1. **Random seeds**: Set for MCMC initialization
2. **Version control**: All package versions recorded
3. **Data versioning**: Checksums for raw data files
4. **Code documentation**: Inline comments and docstrings

## Quality Assurance

### Data Quality Checks
1. **Completeness**: <1% missing climate data
2. **Outlier detection**: ±4 SD flagged for review
3. **Consistency**: Cross-validation with independent sources
4. **Temporal alignment**: Ensure proper date matching

### Model Validation
1. **Posterior predictive checks**: Compare observed vs predicted
2. **Cross-validation**: Leave-one-LHD-out assessment
3. **Residual analysis**: Check for spatial/temporal patterns
4. **Sensitivity to outliers**: Influence diagnostics

## Ethical Considerations

- **Ethics approval**: Not required (aggregate ecological data)
- **Privacy**: No individual-level data used
- **Data governance**: Compliance with NSW Health policies
- **Reporting**: STROBE guidelines for observational studies

## Limitations and Assumptions

### Ecological Fallacy
- Associations at LHD level may not reflect individual risk
- Mitigation: Careful interpretation, acknowledge in discussion

### Exposure Misclassification
- Climate measured at LHD centroid
- Assumes uniform exposure within LHD
- Impact: Likely attenuation of associations

### Unmeasured Confounding
- Socioeconomic factors not explicitly modeled
- Behavioral changes over time
- Addressed partially by model design

### Model Assumptions
1. **Poisson distribution**: Appropriate for count data
2. **Log-linear relationships**: After transformation
3. **Spatial stationarity**: Constant spatial correlation
4. **Missing at random**: For incomplete data

## References

Key methodological references:
1. Gasparrini et al. (2010) - DLNM methodology
2. Besag et al. (1991) - CAR models
3. Rue & Held (2005) - Gaussian Markov random fields
4. de Valpine et al. (2017) - NIMBLE framework