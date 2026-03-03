import os

import glob
import math
import pandas as pd
import numpy as np
from os.path import expanduser 

# =============================
# User settings
# =============================

# Folder containing the AQS CSVs you already generated
CSV_DIR = OUTDIR_CSV = expanduser("~/Documents/Research/AQS_CSV/")   # <-- change if needed

# Site coordinates for filtering
LAT0 = 36.21
LON0 = -81.67
MAX_DIST_KM = 350.0
YEAR = 2024

# =============================
# Haversine distance function
# =============================

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points on Earth in km."""
    R = 6371.0  # Earth radius in km
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (math.sin(dphi/2)**2 +
         math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

    
# =============================
# Step 1: Load all CSVs
# =============================

def combine_and_filter():
    csv_files = glob.glob(os.path.join(CSV_DIR, "*.csv"))
    if not csv_files:
        print("❌ No CSV files found in aqs_csv folder.")
        return

    dfs = []
    for f in csv_files:
        try:
            df = pd.read_csv(f)
            df = read_in_data(df)
            print(f"Loaded {len(df)} rows from {os.path.basename(f)}")
            dfs.append(df)
        except Exception as e:
            print(f"⚠️ Error reading {f}: {e}")

    # Combine
    combo = large_combine(dfs)
    print(f"\nTotal combined rows after filtering: {len(combo)}")


    # Save final CSV
    outname = f"AQS_combined_within_{int(MAX_DIST_KM)}km_{YEAR}.csv"
    outfile = os.path.join(CSV_DIR, outname)    
    combo = combo.drop_duplicates() 
    # filtered = filtered.reset_index()
    # filtered = filtered.set_index(['Datetime(UTC)'])
    combo.to_csv(outfile)

    print(f"\n✅ Saved filtered combined file:")
    print(f"   {outfile}")

def read_in_data(data):
    data['Datetime(UTC)'] = pd.to_datetime(data['date_local'])# + ' ' + data['time_gmt'])
    data['Position'] = np.array(list(map(str, data["latitude"].to_numpy()))) +',' + np.array(list(map(str, data["longitude"].to_numpy())))
    data['Location'] = np.array(list(map(str,data['county'])))+' county,'+np.array(list(map(str,data['state'])))
    units = {
        "Micrograms/cubic meter (LC)": 'µg/m^3 ATP',
        "Micrograms/cubic meter (25 C)": 'µg/m^3 STP',
        "Inverse 100 Megameters" : '0.01 Mm-1',
        "Parts per billion": 'ppb',
        "Parts per million": 'ppm',
    }
    label = f'{data['parameter'].to_numpy()[0]} [{units[data['units_of_measure'].to_numpy()[0]]}]'
    data[label] = data['sample_measurement']
    data = data.set_index(['Datetime(UTC)'])
    kept = ['Location',label,'uncertainty','qualifier','sample_duration','Position','longitude', 'latitude']
    cols = [col for col in data.columns.to_numpy() if col not in kept]
    data = data.drop(columns =cols)
    data = data[kept]
    data["dist_km"] = data.apply(lambda r: haversine_km(r.latitude, r.longitude, LAT0, LON0),axis=1    )
    data = data[data["dist_km"] <= MAX_DIST_KM].copy()
    print(f"Rows WITHIN {MAX_DIST_KM} km: {len(data)}")
    input((data))
    return data

def large_combine(dfs):
    data = dfs[0]
    x_base = list(data.columns.to_numpy()[1].split('['))[0]
    for i in range(len(dfs)-1):
        if i == 0:
            new_base = list(dfs[i+1].columns.to_numpy()[1].split('['))[0]
            data = pd.merge(data, dfs[i+1], how = "outer", left_index=True, right_index=True, suffixes=(f' {x_base}',f' {new_base}'))
            input(data)
        else: 
            new_base = list(dfs[i+1].columns.to_numpy()[1].split('['))[0]
            data = pd.merge(data, dfs[i+1], how = "outer", left_index=True, right_index=True, suffixes=(f'',f' {new_base}'))
            input(data)
    input(data)
    return data


if __name__ == "__main__":
    combine_and_filter()

