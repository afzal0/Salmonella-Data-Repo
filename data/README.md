# Data Documentation

This directory contains all data files used in the analysis of climate impacts on Salmonella infections in NSW.

## Directory Structure

- `raw/`: Original, unmodified data files
- `processed/`: Cleaned and prepared data files ready for analysis

## Data Sources

### 1. Epidemiological Data
- **File**: `raw/datafinal.csv`
- **Source**: NSW Ministry of Health Notifiable Conditions Information Management System (NCIMS)
- **Period**: January 1991 - December 2022
- **Description**: Monthly Salmonella case counts by Local Health District

### 2. Climate Data
- **Included in**: `raw/datafinal.csv`
- **Source**: Bureau of Meteorology (BOM)
- **Variables**:
  - Temperature (mean, max, min)
  - Rainfall (monthly total)
  - Humidity (relative humidity)
- **Processing**: Aggregated from weather stations to LHD level

### 3. Flood Data
- **Included in**: `raw/datafinal.csv`
- **Source**: Emergency Management Australia
- **Variable**: Floodwave (0-10 severity scale)
- **Description**: Monthly flood severity index by LHD

### 4. Spatial Data
- **Files**: 
  - `raw/nswlhd.*` - NSW Local Health District boundaries
  - `raw/STE_2016_AUST.*` - Australian state boundaries
  - `raw/LocalHealthDistricts.kml` - LHD boundaries in KML format
- **Source**: Australian Bureau of Statistics
- **Format**: ESRI Shapefile and KML

### 5. Population Data
- **Included in**: `raw/datafinal.csv`
- **Source**: Australian Bureau of Statistics
- **Description**: Annual population estimates by LHD

## Data Dictionary

### Main Dataset (datafinal.csv)

| Column | Type | Description | Range/Values |
|--------|------|-------------|--------------|
| Date | Date | First day of month | 1991-01-01 to 2022-12-31 |
| LHD | Character | Local Health District name | 15 districts |
| Rainfall | Numeric | Monthly rainfall (mm) | 0 - 500+ |
| Tmean | Numeric | Mean temperature (°C) | 5 - 35 |
| Tmax | Numeric | Maximum temperature (°C) | 10 - 45 |
| Tmin | Numeric | Minimum temperature (°C) | -5 - 30 |
| Elevation | Numeric | Average elevation (m) | 0 - 1500 |
| Floodwave | Integer | Flood severity index | 0 - 10 |
| Salmonella | Integer | Monthly case count | 0 - 200+ |
| Population | Integer | LHD population | 50,000 - 1,000,000+ |
| Incidence | Numeric | Cases per 100,000 | 0 - 50+ |

### Local Health Districts

The 15 LHDs included in the analysis:
1. Central Coast
2. Far West
3. Hunter New England
4. Illawarra Shoalhaven
5. Mid North Coast
6. Murrumbidgee
7. Nepean Blue Mountains
8. Northern NSW
9. Northern Sydney
10. South Eastern Sydney
11. South Western Sydney
12. Southern NSW
13. Sydney
14. Western NSW
15. Western Sydney

## Processed Data Files

### 1. data_casecrossover.RData
- **Created by**: `DataPrep_3.R`
- **Contents**:
  - Cross-basis matrices for DLNM
  - Strata definitions (LHD:Year-Month)
  - Centered exposure variables

### 2. data_timeseries.RData
- **Created by**: `DataPrep_3.R`
- **Contents**:
  - Time series formatted data
  - Temporal basis functions
  - Regional coding

### 3. spatial_adjacency.RData
- **Created by**: `DataPrep_3.R`
- **Contents**:
  - Spatial adjacency matrix
  - Neighbor lists
  - Spatial weights

## Data Quality Notes

1. **Missing Data**: 
   - Minimal missing climate data (<1%)
   - Complete case reporting from 1991 onwards

2. **Data Validation**:
   - Cross-referenced with BOM stations
   - Validated against state health reports

3. **Limitations**:
   - Reporting practices may have changed over time
   - Climate data aggregated to LHD level
   - Flood severity index is semi-quantitative

## Usage Guidelines

1. Always use processed data files for analysis
2. Raw data should not be modified
3. Document any additional data processing steps
4. Check for updates to population estimates

## Contact

For data access questions or issues, contact:
- Data curator: [Name]
- Email: [email address]