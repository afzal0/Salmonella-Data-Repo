# Too Hot, Too Wet: Bayesian Spatial Modelling of Climate-Driven Salmonella Risk in New South Wales, Australia, 1991–2022

This repository contains the data, code, and supplementary materials for replicating the results of our study investigating the relationship between climate variables and Salmonella infections in New South Wales, Australia.

## Authors

Oyelola A. Adegboye¹²*†, Tehan Amarasena²†, Mohammad Afzal Khan³†, Hassan Ajulo², Anton Pak⁴, David Taniar³, Theophilus I. Emeto²

¹ Menzies School of Health Research, Darwin, Charles Darwin University, NT, Australia  
² Public Health and Tropical Medicine, College of Medicine and Dentistry, James Cook University, Townsville, QLD 4811, Australia  
³ Faculty of Information Technology, Monash University, Melbourne, Australia  
⁴ Centre for the Business and Economics of Health, The University of Queensland, Brisbane, Queensland, Australia  

† Equal contribution  
\* Corresponding author: Oyelola Adegboye, Menzies School of Health Research, Darwin, Charles Darwin University, NT, Australia

## Abstract

**Background**: Salmonella infections are responsible for a substantial number of hospitalisations due to bacterial gastrointestinal infections in Australia and contribute significantly to the global burden of infectious diseases. Evidence suggests a seasonal variation of salmonella worldwide, but there is limited literature on these relationships involving rainfall and temperature, particularly in Australia. 

**Objective**: This study investigated the trend of Salmonella cases and their relationship with climatic factors such as temperature, rainfall and flood events at the local health districts (LHDS) level across New South Wales (NSW), Australia. 

**Methods**: Monthly Salmonella cases and climate data (mean temperature, mean minimum temperature, mean maximum temperature and rainfall) were matched to their corresponding LHD in NSW. 

**Results**: Higher recordings of all the climate variables were positively associated with salmonella (all with p < 0.001), and mean monthly minimum temperature was most strongly associated with salmonella (r= 0.404). There is a greater association of mean monthly temperature and total daily rainfall with salmonella in metropolitan LHDs than in the rural/regional LHDs of NSW. 

**Conclusion**: Rising temperatures, increased rainfall, and flood events are likely to contribute to higher rates of salmonella in NSW. These findings have important public health implications that extend beyond NSW, highlighting the need for broader national and global preparedness in the context of a changing climate.

**Keywords**: Salmonella, Foodborne disease, Climate variability, Environmental drivers, Climate change and health, Spatial Bayesian

## Repository Structure

```
.
├── README.md                    # This file
├── data/                        # Data files and documentation
│   ├── raw/                     # Original data files
│   ├── processed/               # Processed data ready for analysis
│   └── README.md               # Data dictionary and sources
├── code/                        # Analysis scripts
│   ├── 01_data_preparation/     # Data preprocessing scripts
│   ├── 02_modeling/             # Statistical modeling scripts
│   ├── 03_visualization/        # Result visualization scripts
│   └── README.md               # Code execution guide
├── results/                     # Generated outputs
│   ├── figures/                 # Plots and maps
│   ├── tables/                  # Summary statistics
│   └── models/                  # Saved model objects
├── docs/                        # Additional documentation
│   ├── requirements.md          # Software dependencies
│   └── methods.md              # Detailed methodology
└── supplementary/               # Supplementary materials

```

## Quick Start

### Prerequisites

- R version 4.2.0 or higher
- RStudio (recommended)
- Required R packages (see `docs/requirements.md`)

### Replication Steps

1. **Clone the repository**
   ```bash
   git clone https://github.com/[username]/salmonella-climate-nsw.git
   cd salmonella-climate-nsw
   ```

2. **Install dependencies**
   ```R
   source("code/00_install_packages.R")
   ```

3. **Run the analysis pipeline**
   ```R
   # Step 1: Data preparation
   source("code/01_data_preparation/DataPrep_3.R")
   
   # Step 2: Model fitting
   source("code/02_modeling/SB_DLNM-New.R")
   
   # Step 3: Generate results
   source("code/03_visualization/DLNM-Curves-2.R")
   source("code/03_visualization/Relative Risk Final.R")
   source("code/03_visualization/Probabibility Maps.R")
   
   # Step 4: Model evaluation
   source("code/02_modeling/Accuracy Assessment.R")
   ```

## Data Description

### Input Data

- **Epidemiological data**: Monthly Salmonella case counts by LHD (1991-2022)
- **Climate data**: Monthly temperature (mean, max, min), rainfall, flood events, and humidity
- **Spatial data**: NSW Local Health District boundaries (shapefile)
- **Population data**: Annual population estimates by LHD

### Key Variables

| Variable | Description | Unit |
|----------|-------------|------|
| Date | Month-year | YYYY-MM-DD |
| LHD | Local Health District name | Text |
| Salmonella | Monthly case count | Count |
| Tmean | Mean temperature | °C |
| Rainfall | Total monthly rainfall | mm |
| Floodwave | Flood severity index | 0-10 scale |
| Humidity | Average relative humidity | % |
| Population | LHD population | Count |
| Incidence | Cases per 100,000 | Rate |

## Statistical Methods

### Models

1. **Model 3: Case-Crossover with Spatial Effects**
   - Conditional Poisson regression
   - Intrinsic CAR spatial correlation
   - Stratification by LHD and year-month

2. **Model 4: Time-Series with Spatial Effects**
   - Poisson regression with temporal trends
   - Seasonal components
   - Regional random intercepts

### Key Features

- Distributed Lag Non-Linear Models (DLNM) to capture delayed effects
- Bayesian inference using NIMBLE
- Spatial correlation using intrinsic Conditional Autoregressive (iCAR) models
- Natural splines for flexible exposure-response relationships

## Results

The analysis produces:

1. **Exposure-response curves**: Non-linear relationships between climate variables and Salmonella risk
2. **Spatial risk maps**: Relative risk by LHD for different exposure levels
3. **Probability maps**: Exceedance probabilities for risk thresholds
4. **Model diagnostics**: Accuracy metrics and model comparison statistics

## Citation

If you use this code or data, please cite:

```
Adegboye, O.A., Amarasena, T., Khan, M.A., Ajulo, H., Pak, A., Taniar, D., & Emeto, T.I. (2025). 
Too Hot, Too Wet: Bayesian Spatial Modelling of Climate-Driven Salmonella Risk in New South Wales, 
Australia, 1991–2022. [Journal Name - in review].
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or issues, please:
- Open an issue on GitHub
- Contact the corresponding author: Oyelola Adegboye at Menzies School of Health Research

## Acknowledgments

We acknowledge the NSW Ministry of Health for providing the epidemiological data and the Bureau of Meteorology for climate data access.