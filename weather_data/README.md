# Weather Data Directory

This directory should contain the climate data files for all NSW Local Health Districts (LHDs) from 1990-2022.

## Required Files

### Primary Data File
- `final_datafile.csv` - **Main dataset** containing:
  - Daily Salmonella case counts by LHD
  - Temperature (mean, max, min)
  - Rainfall measurements
  - Humidity data  
  - Flood wave indicators
  - Population data
  - Date and location identifiers

### Individual LHD Climate Files (Optional)
Individual CSV files for each NSW LHD containing daily climate data:
- `Albury_Wodonga_Health_(Network_with_Victoria).csv`
- `Central_Coast.csv`
- `Far_West.csv`
- `Hunter_New_England.csv`
- `Illawarra_Shoalhaven.csv`
- `Mid_North_Coast.csv`
- `Murrumbidgee.csv`
- `Nepean_Blue_Mountains.csv`
- `Northern_NSW.csv`
- `Northern_Sydney.csv`
- `South_Eastern_Sydney.csv`
- `South_Western_Sydney.csv`
- `Southern_NSW.csv`
- `Sydney.csv`
- `Western_NSW.csv`
- `Western_Sydney.csv`

### Data Processing Scripts (Optional)
- `data_cleaning.py` - Python script for cleaning raw weather data
- `data_wrangle_climate.py` - Climate data wrangling and processing
- `join_data.py` - Script to merge epidemiological and climate data
- `merged_data.csv` - Intermediate merged dataset

## Data Sources

### Climate Data
- **Source**: Australian Bureau of Meteorology
- **Variables**: Temperature, rainfall, humidity
- **Temporal Coverage**: 1990-2022
- **Spatial Coverage**: NSW Local Health Districts

### Epidemiological Data  
- **Source**: NSW Health Department
- **Disease**: Salmonella infections
- **Unit**: Daily case counts by LHD
- **Temporal Coverage**: 1990-2022

### Population Data
- **Source**: Australian Bureau of Statistics
- **Unit**: Annual population estimates by LHD
- **Use**: Calculating incidence rates

## Data Format Requirements

The main `final_datafile.csv` should contain these essential columns:
- `Date` (YYYY-MM-DD format)
- `LHD` (Local Health District name - must match spatial data exactly)
- `Salmonella` (daily case count - laboratory-confirmed infections)
- `Tmean` (mean daily temperature °C)
- `Rainfall` (daily precipitation mm)
- `Humidity` (relative humidity %)
- `Floodwave` (flood indicator variable)
- `Population` (population denominator for incidence calculations)

**Note**: During preprocessing, the analysis creates centered versions:
- `Tmean_centered` = `Tmean` - population mean
- `Rainfall_centered` = `Rainfall` - population mean  
- `Floodwave_centered` = `Floodwave` - population mean
- `Humidity_centered` = `Humidity` - population mean

## Notes for Reproducibility

1. Ensure date formats are consistent (YYYY-MM-DD)
2. Handle missing values appropriately
3. Verify LHD names match spatial data
4. Check for temporal completeness (daily data without gaps)
5. Validate data ranges (temperature, rainfall within realistic bounds)