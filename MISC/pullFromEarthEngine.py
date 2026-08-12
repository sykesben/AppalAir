import ee
import pandas as pd
from datetime import datetime, timezone
from pathlib import Path
from IPython.display import display

# Initialize Earth Engine(make sure its already authenticated)
ee.Initialize(project='appalair-site')
path_out = r'C:\Users\bensy\Documents\Research\EarthEngine'

# Set up a radius identical to EPA AQS dataset
lat = 36.21
lon = -81.69
radius_m = 350000 
# Set time range desired
start_date = '2024-06-01'
end_date = datetime.now(timezone.utc).strftime('%Y-%m-%d')
file_name = fr"\MODIS_weekly_NDVI_FPAR_{start_date}-{end_date}_{radius_m/1000}km.csv"
file_out = Path(path_out+file_name)

# Generate the region to pull from
point = ee.Geometry.Point([lon, lat])
region = point.buffer(radius_m)

# pull the Normalized Difference Vegetation Index
ndvi_ic = ee.ImageCollection('MODIS/061/MOD13A2').filterDate(start_date, end_date).select('NDVI')

# pull the Fraction of Photosynthetically Active Radiation
fpar_ic = ee.ImageCollection('MODIS/061/MCD15A3H').filterDate(start_date, end_date).select('Fpar')

# pull the Fire radiative power
fire_ic = ee.ImageCollection("NASA/VIIRS/002/VNP14A1").filterDate(start_date, end_date).select('MaxFRP')

# pulled burned area in the region -> slightly different process for burn area
burn_ic = ee.ImageCollection("MODIS/061/MCD64A1").filterDate(start_date, end_date).filterBounds(region)

# pull the drought index
drought_ic = ee.ImageCollection("projects/sat-io/open-datasets/us-drought-monitor").select('DM')
# # pull the cloud cover 
# cc_ic = ee.ImageCollection("NOAA/NWS/RTMA").filterDate(start_date, end_date).select('TCDC')

# pull smoke info
smoke_ic = ee.ImageCollection("NOAA/VIIRS/AOD_EDR/V3").select()

#for burn area 
def monthly_metrics(img):
    # Burned area
    burned = img.select('BurnDate').gt(0)
    burned_area = burned.multiply(ee.Image.pixelArea()).reduceRegion(
                reducer=ee.Reducer.sum(),geometry=region,
                scale=500,maxPixels=1e13)
    return ee.Feature(None, {"Date": img.date().format("YYYY-MM"),"burned[km2]":ee.Number(burned_area.get("BurnDate")).divide(1e6)})

features = burn_ic.map(monthly_metrics)
feature_list = features.getInfo()['features']
records = []
for f in feature_list:
    records.append(f['properties'])
df = pd.DataFrame(records)
df['Date'] = pd.to_datetime(df['Date'])
df = df.set_index('Date')

print('data selected')
# Generate monthly formats
months = pd.date_range(start=start_date,end=end_date,freq='MS')
weeks = pd.date_range(start=start_date,end=end_date,freq='W')
rows = []
for start in weeks:
    print(f'pulling for {start}')
    end = start + pd.offsets.DateOffset(weeks=1)
    start_str = start.strftime('%Y-%m-%d')
    end_str = end.strftime('%Y-%m-%d')
    # Monthly composites
    ndvi_img = ndvi_ic.filterDate(start_str, end_str).mean()
    fpar_img = fpar_ic.filterDate(start_str, end_str).mean()
    # burn_img = burn_ic.filterDate(start_str, end_str).mean()
    fire_img = fire_ic.filterDate(start_str, end_str).mean()
    drought_img = drought_ic.filterDate(start_str, end_str).mean()
    # cc_img = cc_ic.filterDate(start_str, end_str).mean()
    # Regional averages
    ndvi_mean = ndvi_img.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=1000,
        maxPixels=1e13)
    fpar_mean = fpar_img.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=500,
        maxPixels=1e13)
    # burn_mean = burn_img.reduceRegion(
    #     reducer=ee.Reducer.sum(),
    #     geometry=region,
    #     scale=1000,
    #     maxPixels=1e13)
    fire_mean = fire_img.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=region,
        scale=1000,
        maxPixels=1e13)
    drought_mean = drought_img.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=250,
        maxPixels=1e13)
    # cc_mean = cc_img.reduceRegion(
    #     reducer=ee.Reducer.mean(),
    #     geometry=region,
    #     scale=500,
    #     maxPixels=1e13)
    try:
        ndvi_val = ndvi_mean.get('NDVI').getInfo()
        ndvi_val *= 0.0001
    except:
        ndvi_val = None
    try:
        fpar_val = fpar_mean.get('Fpar').getInfo()
        fpar_val *= 0.01
    except:
        fpar_val = None
    # try:
    #     burn_val = burn_mean.get('MaxFRP').getInfo()
    #     burn_val *= 0.1
    # except: 
    #     burn_val = None
    try:
        fire_val = fire_mean.get('MaxFRP').getInfo()
    except: 
        fire_val = None
    try: 
        drought_val = drought_mean.get('DM').getInfo()
    except:
        drought_val = None
        # try: 
    #     cc_val = cc_mean.get('TCDC').getInfo()
    # except:
    #     cc_val = None

    rows.append({
        'Date': start.strftime('%Y-%m-%d'),
        'NDVI': ndvi_val,
        'FPAR': fpar_val,
        # 'MaxFRP_MODUS': burn_val,
        'MaxFRP[MW]': fire_val,
        'Drought' : drought_val
        })

# Save to CSV
df_tot = pd.DataFrame(rows)
df_tot['Date'] = pd.to_datetime(df_tot['Date'])
df_tot = df_tot.set_index('Date')
# df_tot = pd.merge(df_tot,df,left_index=True,right_index=True)
print(df_tot.head())
df_tot.to_csv(file_out)

