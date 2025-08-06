import pandas as pd

# Read the main data file (datafinal file)
datafinal = pd.read_csv("datafinal.csv", parse_dates=["Date"])

# Read the merged humidity data file (with LHD, Date, Humidity)
humidity = pd.read_csv("merged_data.csv", parse_dates=["Date"])

# Merge the two datasets on LHD and Date using a left join
merged_data = pd.merge(datafinal, humidity, on=["LHD", "Date"], how="left")

# Debug: Count the total records and how many records did not join (missing Humidity)
total_datafinal = len(datafinal)
total_humidity = len(humidity)
total_merged = len(merged_data)
missing_humidity_count = merged_data['Humidity'].isna().sum()

print(f"Total records in datafinal file: {total_datafinal}")
print(f"Total records in humidity file: {total_humidity}")
print(f"Total records after merge: {total_merged}")
print(f"Records with missing Humidity after merge: {missing_humidity_count}")

# If there are missing records, check which LHD values are affected
if missing_humidity_count > 0:
    missing_lhd = merged_data.loc[merged_data['Humidity'].isna(), 'LHD'].unique()
    print("LHD(s) with missing humidity data:", missing_lhd)

# Save the merged DataFrame to a new CSV file
merged_data.to_csv("final_datafile.csv", index=False)
print("Merged data saved as final_datafile.csv")
