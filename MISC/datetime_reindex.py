import os

import glob
import math
import pandas as pd
import numpy as np
from os.path import expanduser 


def read_in_data(data):
    data['Datetime(UTC)'] = pd.to_datetime(data['date_gmt'] + ' ' + data['time_gmt'])
    data['Position'] = np.array(list(map(str, data["latitude"].to_numpy()))) +',' + np.array(list(map(str, data["longitude"].to_numpy())))
    data['Location'] = np.array(list(map(str,data['county'])))+' county,'+np.array(list(map(str,data['state'])))
    units = {
        "Micrograms/cubic meter (LC)": '[µg/m^3 ATP]',
        "Micrograms/cubic meter (25 C)": '[µg/m^3 STP]',
        "Inverse 100 Megameters" : '[0.01 Mm-1]',
        "Parts per billion": '[ppb]',
        "Parts per million": '[ppm]',
    }
    label = [f"{a} {b}" for a, b in zip(data['parameter'].to_numpy(),np.array([units[key] for key in data['units_of_measure'].to_numpy() if key in units]))]
    data['label'] = label
    data = data.set_index(['Datetime(UTC)','Position'])
    kept = ['Location','label','sample_measurement','uncertainty','qualifier','sample_duration','method', 'dist_km'] #observation_percent
    cols = [col for col in data.columns.to_numpy() if col not in kept]
    data = data.drop(columns =cols)
    data = data[kept]
    data = data.sort_index()
    daily = data[data['sample_duration']=='24 HOUR']
    cont = data[data['sample_duration']=='1 HOUR']
    daily = daily.drop(columns='sample_duration')
    cont = cont.drop(columns='sample_duration')
    cont = cont.reset_index()
    cont = pd.pivot_table(cont, values='sample_measurement', index=['Datetime(UTC)'],
                       columns=['Position'], aggfunc="mean")
    cont['Measurement'] ='PM2.5 [µg/m^3 ATP]'
    daily = daily.reset_index()
    daily['Date(UTC)'] = daily['Datetime(UTC)'].dt.date
    daily = pd.pivot_table(daily, values='sample_measurement', index=['Date(UTC)','Position','Location',],
                       columns=['label'], aggfunc="mean")
    input(daily)
    data.to_csv(r"C:\Users\bensy\Documents\Research\AQS_CSV\AQS_combined_2025.csv")
    daily.to_csv(r"C:\Users\bensy\Documents\Research\AQS_CSV\AQS_combined_daily_2025.csv")
    cont.to_csv(r"C:\Users\bensy\Documents\Research\AQS_CSV\AQS_combined_cont_2025.csv")
    return data

f = r"C:\Users\bensy\Documents\Research\AQS_CSV\AQS_combined_within_350km_2025.csv"
big_file = pd.read_csv(f)


if __name__ == "__main__":
    read_in_data(big_file)
