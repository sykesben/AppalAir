from activationFraction import smps_data, cpc_data, ccn_data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

cpc = r"C:\Users\bensy\Documents\Research\CPC_2026.csv"
smps = [r"C:\Users\bensy\Documents\Research\SMPS_2026\2026_SMPS_NumberSizeDist_1hr_clean.csv"]
CCN = [r"C:\Users\bensy\Documents\Research\CCN\Processed\CCN_lvl2_2026_1hr.csv"]

smpsDf, cols = smps_data(smps, 'h')
print(cols)
cpcDF = cpc_data(cpc,'h')
ccndf,ccn_cols,ss_cols = ccn_data(CCN,'h')
df =pd.merge(cpcDF,smpsDf[cols],left_index = True, right_index = True, how='right')
df =pd.merge(df,ccndf,left_index = True, right_index = True)
df = df.dropna(axis='index', how ='all')
df['Fact CPC'] = df['N(cm-3)_avg_setpt0.7']/df['CPC N[cm-3]']
df['Fact SMPS'] = df['N(cm-3)_avg_setpt0.7']/df['Total Concentration (#/cm³)']
plt.ion()
fig, ax = plt.subplots()
ax.plot(df['Fact CPC'],label ='Fact CPC')
ax.plot(df['Fact SMPS'] ,label ='Fact SMPS')
ax.legend()
ax.set_ylabel('F[1]')
ax.set_xlabel('Date')
ax.set_title('Activation Fraction Comparisons')
input('Press enter to exit plot...')
plt.ioff()