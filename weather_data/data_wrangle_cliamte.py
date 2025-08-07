import openmeteo_requests
import requests_cache
import pandas as pd
from retry_requests import retry

# Setup the Open-Meteo API client with cache and retry on error
cache_session = requests_cache.CachedSession('.cache', expire_after=-1)
retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
openmeteo = openmeteo_requests.Client(session=retry_session)

# List of locations
locations = [
    # {"location_id": 5, "latitude": -33.399998, "longitude": 151.4, "elevation": 44, "LHD": "Central Coast"},
    # {"location_id": 6, "latitude": -31.8, "longitude": 142.80002, "elevation": 178, "LHD": "Far West"},
    # {"location_id": 7, "latitude": -30.7, "longitude": 150.80002, "elevation": 415, "LHD": "Hunter New England"},
    # {"location_id": 8, "latitude": -35, "longitude": 150.4, "elevation": 262, "LHD": "Illawarra Shoalhaven"},
    # {"location_id": 9, "latitude": -30.899998, "longitude": 152.70001, "elevation": 94, "LHD": "Mid North Coast"},
    # {"location_id": 10, "latitude": -34.7, "longitude": 146.20001, "elevation": 140, "LHD": "Murrumbidgee"},
    # {"location_id": 0, "latitude": -33.399998, "longitude": 150.4, "elevation": 728, "LHD": "Nepean Blue Mountains"},
    # {"location_id": 12, "latitude": -29.2, "longitude": 152.9, "elevation": 91, "LHD": "Northern NSW"},
    # {"location_id": 1, "latitude": -33.7, "longitude": 151.1, "elevation": 121, "LHD": "Northern Sydney"},
    # {"location_id": 15, "latitude": -33.95402259, "longitude": 151.3425526, "elevation": 6, "LHD": "South Eastern Sydney"},
    # {"location_id": 2, "latitude": -34, "longitude": 151.5, "elevation": 0, "LHD": "South Western Sydney"},
    # {"location_id": 13, "latitude": -35.8, "longitude": 149.4, "elevation": 1125, "LHD": "Southern NSW"},
    # {"location_id": 3, "latitude": -34.2, "longitude": 150.5, "elevation": 416, "LHD": "Sydney"},
    {"location_id": 14, "latitude": -31.399998, "longitude": 147.20001, "elevation": 166, "LHD": "Western NSW"},
    {"location_id": 4, "latitude": -33.899998, "longitude": 151.1, "elevation": 28, "LHD": "Western Sydney"}
]

# Function to get weather data for a location
def get_weather_data(latitude, longitude):
    url = "https://archive-api.open-meteo.com/v1/archive"
    params = {
        "latitude": latitude,
        "longitude": longitude,
        "start_date": "1990-01-01",
        "end_date": "2022-12-31",
        "hourly": ["relative_humidity_2m"],
        "timezone": "auto"
    }
    return openmeteo.weather_api(url, params=params)

# Loop over each location, download data, and save to CSV
for loc in locations:
    response = get_weather_data(loc['latitude'], loc['longitude'])[0]
    
    hourly = response.Hourly()
    hourly_data = {
        "date": pd.date_range(
            start=pd.to_datetime(hourly.Time(), unit="s", utc=True),
            end=pd.to_datetime(hourly.TimeEnd(), unit="s", utc=True),
            freq=pd.Timedelta(seconds=hourly.Interval()),
            inclusive="left"
        ),
        "relative_humidity_2m": hourly.Variables(0).ValuesAsNumpy(),
    }
    
    # Create DataFrame
    df = pd.DataFrame(hourly_data)
    
    # Save to CSV using the LHD name
    filename = f"{loc['LHD'].replace(' ', '_')}.csv"
    df.to_csv(filename, index=False)

# Confirmation of completion
print("Weather data for all locations downloaded and saved to CSV.")
