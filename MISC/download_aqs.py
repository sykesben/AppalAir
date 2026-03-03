import os
import time
import json
from os.path import expanduser 

import requests
import pandas as pd

# =============================
# User settings
# =============================

EMAIL = "sykesbb@appstate.edu"  # <-- put your AQS-registered email here
KEY = "saffroncrane32"         # <-- put your AQS key here

BASE_URL = "https://aqs.epa.gov/data/api/sampleData/byState"
BDATE = "20250101"
EDATE = "20250131"


# State codes and species (parameters)
states = {
    "NC": "37"}
'''
    "GA": "13",
    "KY": "21",
    "SC": "45",
    "TN": "47",
    "VA": "51",
    "WV": "54",
'''

params = {
    "NH4":  "88301",
    "EC":   "88321",
    "NO3":  "88306",
    "OC":   "88320",
    "SO4":  "88403",
    "O3": '44201',
    "NOx" : "42600",
    "Light Absorption Coeff":"63102"}
    # "PM25": "88101",  
    # "PM10": "81102",

# Output directories
OUTDIR_JSON = expanduser("~/Documents/Research/AQS_Json/")
OUTDIR_CSV = expanduser("~/Documents/Research/AQS_CSV/")
os.makedirs(OUTDIR_JSON, exist_ok=True)
os.makedirs(OUTDIR_CSV, exist_ok=True)

# Optional: small delay between API calls (seconds)
SLEEP_SEC = 1


def fetch_aqs_response(state_code, param_code):
    """Fetch full JSON response from AQS sampleData/byState for one state & parameter."""
    params_req = {
        "email": EMAIL,
        "key": KEY,
        "param": param_code,
        "bdate": BDATE,
        "edate": EDATE,
        "state": state_code,
       # "format": "application/json",
    }

    resp = requests.get(BASE_URL, params=params_req, timeout=300)

    try:
        resp.raise_for_status()
    except requests.HTTPError as e:
        print(f"HTTP error for state={state_code}, param={param_code}: {e}")
        return None

    try:
        data = resp.json()
    except json.JSONDecodeError as e:
        print(f"JSON decode error for state={state_code}, param={param_code}: {e}")
        return None

    # Print header status if present
    if "Header" in data and data["Header"]:
        status = data["Header"][0].get("status", "unknown")
        msg = data["Header"][0].get("message", "")
        print(f"API status: {status} - {msg}")

    return data


def main():
    for state_abbr, state_code in states.items():
        for species_name, param_code in params.items():
            print(f"\nQuerying: state={state_abbr} ({state_code}), parameter={species_name} ({param_code})")

            data = fetch_aqs_response(state_code, param_code)
            if data is None:
                print("  Skipping (no response or error).")
                time.sleep(SLEEP_SEC)
                continue

            # ----------------------------
            # 1) Save full JSON response
            # ----------------------------
            base_name = f"AQS_{species_name}_{param_code}_{state_abbr}_{state_code}_{BDATE}_{EDATE}"
            json_path = os.path.join(OUTDIR_JSON, base_name + ".json")

            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)

            print(f"  Saved JSON to {json_path}")

            # ------------------------------------------------------
            # 2) Convert Data portion to CSV (like Excel would after
            #    loading JSON and expanding the 'Data' list)
            # ------------------------------------------------------
            if "Data" not in data or not data["Data"]:
                print("  No 'Data' records in JSON; CSV not created.")
            else:
                # data["Data"] is already a list of record dicts
                df = pd.DataFrame(data["Data"])

                csv_path = os.path.join(OUTDIR_CSV, base_name + ".csv")
                df.to_csv(csv_path, index=False)
                print(f"  Saved {len(df)} data rows to {csv_path}")

            time.sleep(SLEEP_SEC)


if __name__ == "__main__":
    main()

