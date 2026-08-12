import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path, PurePath
from datetime import datetime

def smps_means(vol,mass,freq='d'):
    '''
    Takes in a list of smps files and returns some values
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    smps = pd.DataFrame()
    for i in range(len(vol)): #read in smps files and combine
        m = mass[i]
        v = vol[i]
        mfile =pd.read_csv(m) #read in smps file
        vfile= pd.read_csv(v)
        try: 
            mfile = mfile.set_index('Datetime(UTC)')
            vfile = vfile.set_index('Datetime(UTC)')
        except:
            mfile = mfile.set_index("DateTime Sample Start") #Set index
            vfile = vfile.set_index("DateTime Sample Start") 
        mdata = mfile['Total Concentration (µg/m³)'].to_frame()
        vdata = vfile['Total Concentration (nm³/cm³)'].to_frame()
        data = pd.merge(mdata,vdata,right_index=True,left_index=True)
        data['Density [g/cm³]'] = data['Total Concentration (µg/m³)']/data['Total Concentration (nm³/cm³)']*(10**9)
        smps = data if smps.empty else pd.concat([smps, data])
    smps.index = pd.to_datetime(smps.index, format='mixed')
    smps = smps.resample(freq).mean()
    return smps

mass = [r"C:\Users\bensy\Documents\Research\SMPS_2026\2026_SMPS_MassSizeDist_1hr_clean_stp.csv"]
vol = [r"C:\Users\bensy\Documents\Research\SMPS_2026\2026_SMPS_VolumeSizeDist_1hr_clean_stp.csv"]
smps = smps_means(vol,mass)
input(smps)