

import os
import glob
import math
import pandas as pd

# =============================
# User settings
# =============================

# Folder containing the AQS CSVs you already generated
CSV_DIR = r'C:\Users\bensy\Documents\Research\AQS_CSV'

# Site coordinates for filtering
LAT0 = 36.21
LON0 = -81.69
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
            print(f"Loaded {len(df)} rows from {os.path.basename(f)}")
            dfs.append(df)
        except Exception as e:
            print(f"⚠️ Error reading {f}: {e}")

    # Combine
    combined = pd.concat(dfs, ignore_index=True)
    print(f"\nTotal combined rows BEFORE filtering: {len(combined)}")

    # Check required columns
    if "latitude" not in combined.columns or "longitude" not in combined.columns:
        raise ValueError("❌ Expected columns 'latitude' and 'longitude' not found.")

    # Clean coordinate columns
    combined["latitude"] = pd.to_numeric(combined["latitude"], errors="coerce")
    combined["longitude"] = pd.to_numeric(combined["longitude"], errors="coerce")
    combined = combined.dropna(subset=["latitude", "longitude"])

    # Compute distances
    combined["dist_km"] = combined.apply(
        lambda r: haversine_km(r.latitude, r.longitude, LAT0, LON0),
        axis=1
    )

    # Filter
    filtered = combined[combined["dist_km"] <= MAX_DIST_KM].copy()
    print(f"Rows WITHIN {MAX_DIST_KM} km: {len(filtered)}")

    # Save final CSV
    outname = f"AQS_combined_within_{int(MAX_DIST_KM)}km_{YEAR}.csv"
    outfile = os.path.join(CSV_DIR, outname)
    filtered.to_csv(outfile, index=False)

    print(f"\n✅ Saved filtered combined file:")
    print(f"   {outfile}")


if __name__ == "__main__":
    combine_and_filter()

