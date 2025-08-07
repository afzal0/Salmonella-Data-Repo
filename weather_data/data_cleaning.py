import pandas as pd
import glob
import os

# Get a list of all CSV files in the current directory
csv_files = glob.glob("*.csv")

merged_data = []

for file in csv_files:
    # Extract LHD name from file name (remove extension and replace underscores with spaces)
    lhd_name = os.path.splitext(os.path.basename(file))[0].replace("_", " ")
    
    # Read CSV file; the date column in your files is named "date"
    df = pd.read_csv(file, parse_dates=["date"])
    
    # Convert the 'date' column to the first day of its month
    df['MonthStart'] = df['date'].dt.to_period('M').dt.to_timestamp()
    
    # Rename the relative humidity column to "Humidity"
    df.rename(columns={"relative_humidity_2m": "Humidity"}, inplace=True)
    
    # Aggregate hourly data to monthly by calculating the mean humidity
    df_monthly = df.groupby('MonthStart', as_index=False)['Humidity'].mean()
    
    # Add the LHD column using the extracted file name
    df_monthly['LHD'] = lhd_name
    
    # Rename 'MonthStart' to 'Date' and reorder the columns
    df_monthly.rename(columns={'MonthStart': 'Date'}, inplace=True)
    df_monthly = df_monthly[['LHD', 'Date', 'Humidity']]
    
    merged_data.append(df_monthly)

# Combine data from all CSV files into one DataFrame
merged_df = pd.concat(merged_data, ignore_index=True)

# Optional: sort by LHD and Date
merged_df.sort_values(['LHD', 'Date'], inplace=True)

# Save the merged DataFrame to a new CSV file
merged_df.to_csv("merged_data.csv", index=False)

print("Merged CSV file saved as merged_data.csv")
