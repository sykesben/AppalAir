"""
Date: 4/25/2026
Author: Ben Sykes
Purpose: Call comb_files to combine data from different datasets
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import least_squares as LSfit
from datetime import datetime 
from scipy.interpolate import interp1d
import re

def master_data(f,freq='d'):
    '''
    Takes in the master file and specifically cuts out the AQS data
    ----------

    Parameters
    ++++++++++
    f : [list of str] Paths to Master file
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    master : [DataFrame] Master data file
    spec : [list of str] Names of used columns from chemistry output
    '''
    master=pd.read_csv(f) #read in AQS file
    master=master.set_index("Local time (UTC-5)") #Set index
    master['Date(UTC)'] = pd.to_datetime(master.index) + pd.Timedelta(hours=5)
    specs = ['NH4_11000','SO4_11000','NO3_11000','Org_11000','1hrMC_µg/m3','org/total','SO4/total'] #important speciated data columns
    master=master.set_index("Date(UTC)")
    master = master[specs]
    master.columns = master.columns.str.replace('_11000', ' [µg/m3] ACSM')
    master = master.resample(freq).mean()
    master = master.dropna()
    return master,master.columns.to_numpy()

def ACSM_data(f, freq = 'd'):
    acsm = pd.read_excel(f)
    print('0')
    acsm = acsm.set_index('Local time (UTC-5)')
    acsm.index = pd.to_datetime(acsm.index, format='mixed')+pd.Timedelta('5h')
    acsm.index.rename('Datetime(UTC)', inplace=True)
    acsm = acsm[['Chl_11000','NH4_11000','SO4_11000','NO3_11000','Org_11000',
                 'SO4 / Org','NO3 / Org','NH4/ Org','NO3 / SO4','org/total','SO4/total','NO3/total','Total mass',
                 'PPPMF_OOA','PPPMF_HOA','f_org43','f_org44','f_org60']]
    acsm.columns = acsm.columns.str.replace('_11000', '[ug/m3]')
    acsm[acsm<0] = 0
    with pd.option_context('display.max_rows', 5, 'display.max_columns', None):
        input(acsm)

print('0')
f = r"C:\Users\bensy\Documents\Research\20260512_ACSM_240529-260511_speciated_allmz(1).xlsx"
ACSM_data(f)