import pandas as pd
import geopandas as gpd
from datetime import datetime, timezone
from pathlib import Path
from IPython.display import display
from shapely import Point
import matplotlib.pyplot as plt
import numpy as np

hms26 = r"C:\Users\bensy\Documents\Research\SmokeData\hms_smoke2026.shp"
hms25 = r"C:\Users\bensy\Documents\Research\SmokeData\hms_smoke2025.shp"
hms24 = r"C:\Users\bensy\Documents\Research\SmokeData\hms_smoke2024.shp"

lat = 36.21
lon = -81.69
radius_m = 1_000
region = Point(lat,lon).buffer(radius_m)

# Set time range desired
start_date = '2024-06-01'
end_date = datetime.now(timezone.utc).strftime('%Y-%m-%d')

hms24_data = gpd.read_file(hms24)
hms25_data = gpd.read_file(hms25)
hms26_data = gpd.read_file(hms26)

hms =  gpd.GeoDataFrame(pd.concat([hms24_data,hms25_data, hms26_data],ignore_index=True), crs=gpd.read_file(hms24).crs)
hms_flt = hms.loc[hms.within(region)]
hms_flt["Start"] = pd.to_datetime(hms_flt["Start"], format="%Y%j %H%M")
hms_flt["End"] = pd.to_datetime(hms_flt["End"], format="%Y%j %H%M")
hms_mode = hms_flt.groupby('Start')['Density'].agg(lambda x: pd.Series.mode(x)[0]).to_frame()
hms_mode = hms_mode.resample('d').first()
# dens_dict = {'Heavy':3,'Medium':2,'Light':1}
# hms_mode['Density'] = hms_mode['Density'].apply(lambda x : dens_dict[x])
hms_mode.to_csv(r"C:\Users\bensy\Documents\Research\SmokeData\dailySmoke.csv")