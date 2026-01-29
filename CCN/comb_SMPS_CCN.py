"""
Date: 1/20/26
Author: Ben Sykes
Purpose: Read through SMPS and CCN data to combine
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import os
import datetime as dt 
import scipy as sp
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
from pathlib import Path
from os.path import expanduser 

#Read in files, feel free to replace these with exact 
ccn24 = pd.read_csv(expanduser("~/Documents/Research/CCN_Processed_2024_1hr.csv"))
ccn25 = pd.read_csv(expanduser("~/Documents/Research/CCN_Processed_2025_1hr.csv"))
smps24 =pd.read_csv(expanduser("~/Documents/Research/2024NumConcAVG.csv"))
smps25 =pd.read_csv(expanduser("~/Documents/Research/SMPS_NumberSizeDist_2025_1hr.csv"))

smps25=smps25.set_index("DateTime Sample Start")
smps24 = smps24.set_index('DateTime Sample Start')
ccn24=ccn24.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
ccn25=ccn25.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
ccn24.index = pd.to_datetime(ccn24.index)
ccn25.index = pd.to_datetime(ccn25.index, yearfirst=True)

ccn = pd.concat([ccn24,ccn25])
smps =pd.concat([smps24,smps25])
ccn.index = pd.to_datetime(ccn.index) #format='%Y-%m-%d %H:%M:%S'
smps.index = pd.to_datetime(smps.index)

ccn = ccn[['N(cm-3)_avg_setpt0.1', 'N(cm-3)_avg_setpt0.7','T(C)_inlet','T1(C)','T(C)_sample','T(C)_OPC','T(C)_nafion','Q(lpm)_sample','Q(lpm)_sheath','P(hPA)_sample','ss(%)_calc_setpt0.1', 'ss(%)_calc_setpt0.7',]]
# ccn = ccn.resample('d').mean()
ccn.index.names = ['Date']

numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]
gr200 = [s for s in numsmps if float(s) >200]
gr80 = [s for s in numsmps if float(s) >80]
smps = smps[numsmps]
# smps = smps.resample('d').mean()
smps.index.names = ['Date']
smps['>80nm'] = smps[gr80].mean(axis=1)  
smps['>200nm'] =smps[gr200].mean(axis=1)


data = pd.merge(ccn,smps[['>80nm', '>200nm']],left_index = True, right_index = True)

data = data.loc[pd.to_datetime('2025-07-25'):pd.to_datetime('2025-11-01')]
data.to_csv(expanduser("~/Documents/Research/SMPS_CCN_comparison_hourly.csv"))
print('done')