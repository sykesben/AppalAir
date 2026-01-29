"""
Date: 10/31/2025
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

ccn17 = pd.read_csv(r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr_17.csv")
ccn16 = pd.read_csv(r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr_16.csv")
ccnfit = pd.read_csv(r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr_fit.csv")
smps =pd.read_csv(r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv")

ccn17=ccn17.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
ccn16=ccn16.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
ccnfit=ccnfit.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
smps=smps.set_index("DateTime Sample Start")
ccn17.index = pd.to_datetime(ccn17.index)
ccn16.index = pd.to_datetime(ccn16.index)
ccnfit.index = pd.to_datetime(ccnfit.index)
smps.index = pd.to_datetime(smps.index)

ccn17 = ccn17[['N(cm-3)_cor_setpt0.1', 'N(cm-3)_cor_setpt0.7']]
ccn17 = ccn17.resample('d').mean()
ccn17.index.names = ['Date']

ccn16 = ccn16[['N(cm-3)_cor_setpt0.1', 'N(cm-3)_cor_setpt0.7']]
ccn16 = ccn16.resample('d').mean()
ccn16.index.names = ['Date']

ccnfit = ccnfit[['N(cm-3)_cor_setpt0.1', 'N(cm-3)_cor_setpt0.7']]
ccnfit = ccnfit.resample('d').mean()
ccnfit = ccnfit.add_suffix('_fit', axis='columns')
ccnfit.index.names = ['Date']

numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]
gr200 = [s for s in numsmps if float(s) >200]
gr100 = [s for s in numsmps if float(s) >100]
gr90 = [s for s in numsmps if float(s) >90]
gr80 = [s for s in numsmps if float(s) >80]
smps = smps[numsmps]
smps = smps.resample('d').mean()
smps.index.names = ['Date']
smps['>80nm'] = smps[gr80].mean(axis=1)  
smps['>90nm'] = smps[gr90].mean(axis=1)  
smps['>100nm'] = smps[gr100].mean(axis=1)  
smps['>200nm'] =smps[gr200].mean(axis=1)


ccn = pd.merge(ccn16, ccn17, left_index=True, right_index=True, suffixes=['_16','_17'])
ccn = pd.merge(ccn,ccnfit, left_index=True, right_index=True)
data = pd.merge(ccn,smps[['>80nm','>90nm','>100nm','>200nm']],left_index = True, right_index = True)
data.index = pd.to_datetime(data.index)
data = data.loc[(data != 0).any(axis=1)]
data = data.dropna()
data= data.resample('ME').mean()
data.to_csv(r"C:\Users\bensy\Documents\Research\SMPS_CCN_comparison monthly.csv")
print('Done')