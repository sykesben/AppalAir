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
from AQS_ACSM_comb import PM25_data
import re
from plotgen import monthly_box_call, hourly_box_call,line_plot, scat_call

def master_data(f,freq='h'):
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

def smps_vol(files,freq='d'):
    '''
    Takes in a list of smps volume files and return total volume
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] 
    '''
    smps = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        try: 
            file = file.set_index('Datetime(UTC)')
        except:
            file=file.set_index("DateTime Sample Start") #Set index
        smps = file if smps.empty else pd.concat([smps, file])
    smps.index = pd.to_datetime(smps.index, format='mixed')
    numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]
    nums = smps[numsmps].copy()
    Dp = nums.columns.astype(float).to_numpy()
    smps = smps[smps['Total Concentration (nm³/cm³)'].notna()]
    smps = smps['Total Concentration (nm³/cm³)']
    return smps

def AE33_data(files,freq = 'h'):
    ae33 = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        file=file.set_index("DateTimeUTC") #Set index
        ae33 = file if ae33.empty else pd.concat([ae33, file])
    ae33.index = pd.to_datetime(ae33.index, format='mixed')
    ae33.index.rename('Datetime(UTC)', inplace=True)
    ae33 = ae33.resample(freq).mean()
    # input(ae33)
    return ae33

def ACSM_data(f, freq = 'h',BC ='',vol = ''):
    MW = {'SO4':96.0626*10**6, 'NH4':18.03846*10**6, 'NO3':62.0049*10**6, 'Chl':35.4530*10**6} #[ug/mol]
    acsm = pd.read_excel(f)
    acsm = acsm.set_index('Local time (UTC-5)')
    acsm.index = pd.to_datetime(acsm.index, format='mixed')+pd.Timedelta('5h')
    acsm.index.rename('Datetime(UTC)', inplace=True)
    print(freq)
    acsm = acsm.resample(freq).mean()
    acsm = acsm[['Chl_11000','NH4_11000','SO4_11000','NO3_11000','Org_11000',
                 'SO4 / Org','NO3 / Org','NH4/ Org','NO3 / SO4','org/total','SO4/total','NO3/total','Total mass',
                 'PPPMF_OOA','PPPMF_HOA','f_org43','f_org44','f_org60']]
    acsm.columns = acsm.columns.str.replace('_11000', '[ug/m3]')
    acsm.columns = acsm.columns.str.replace('/ ','/')
    acsm.columns = acsm.columns.str.replace(' /','/')

    acsm['OOA/OA'] = acsm['PPPMF_OOA']/acsm['Org[ug/m3]']
    acsm['HOA/OA'] = acsm['PPPMF_HOA']/acsm['Org[ug/m3]']
    acsm['OOA/OA'][acsm['OOA/OA']>1] = 1
    acsm['HOA/OA'][acsm['HOA/OA']>1] = 1
    acsm['PM2.5 [ug/m3]'] = acsm['Total mass']
    acsm['O/C'] = 3.82*acsm['f_org44'] + 0.0794
    acsm['H/C'] = 1.01 + 6.07*acsm['f_org43'] - (acsm['f_org44']**2)*16.01
    if (isinstance(BC, pd.Series)):
        acsm = pd.merge(acsm,BC, how='left', left_index=True,right_index=True)
        # acsm['BC[ug/m3]'] = acsm['BC[ug/m3]'].ffill().bfill()
    acsm['k_org'] = 0.18*acsm['O/C']+0.03
    acsm['rho_org'] = (12+acsm['H/C'].fillna(0)+16*acsm['O/C'])/(7+5*acsm['H/C'].fillna(0)+4.15*acsm['O/C']) #[g/cm³]
    acsm['k_org'] = acsm['k_org'].fillna(0.1)#[acsm['k_org']<0.5]
    acsm['rho_org'] = acsm['rho_org'].fillna(1.0)
    acsm['V_org [nm3/cm3]'] = acsm['Org[ug/m3]']/(acsm['rho_org'])*(10**9)
    acsm['V_SO4 [nm3/cm3]'] = acsm['SO4[ug/m3]']/(1.77)*(10**9) 
    acsm['V_NO3 [nm3/cm3]'] = acsm['NO3[ug/m3]']/(1.72)*(10**9)
    acsm['V_NH4 [nm3/cm3]'] = acsm['NH4[ug/m3]']/(1.74)*(10**9)
    acsm = acsm.dropna(axis='rows',how='all')
    # acsm.fillna(value=0,inplace=True)
    if (isinstance(BC, pd.Series)):
        acsm['V_BC [nm3/cm3]'] = acsm['BC[ug/m3]']/(1.8)*(10**9)
        input(acsm[['V_org [nm3/cm3]','V_BC [nm3/cm3]']])

        acsm['V_tot [nm3/cm3]'] = (acsm['V_NH4 [nm3/cm3]'].fillna(0)+acsm['V_SO4 [nm3/cm3]'].fillna(0)+acsm['V_NO3 [nm3/cm3]'].fillna(0)
                                   +acsm['V_org [nm3/cm3]'].fillna(0)+acsm['V_BC [nm3/cm3]'].fillna(0))
        acsm['V_test [nm3/cm3]'] = (acsm['V_NH4 [nm3/cm3]'].fillna(0)+acsm['V_SO4 [nm3/cm3]'].fillna(0)+acsm['V_NO3 [nm3/cm3]'].fillna(0)
                                    +acsm['V_org [nm3/cm3]'].fillna(0))
        # input(acsm[['V_test [nm3/cm3]','V_tot [nm3/cm3]']])
        acsm['MAF'] = 1- acsm['V_BC [nm3/cm3]']/acsm['V_tot [nm3/cm3]']
    else:
        acsm['V_tot [nm3/cm3]'] = (acsm['V_NH4 [nm3/cm3]'].fillna(0)+acsm['V_SO4 [nm3/cm3]'].fillna(0)+acsm['V_NO3 [nm3/cm3]'].fillna(0)
                                    +acsm['V_org [nm3/cm3]'].fillna(0))
    if (isinstance(vol, pd.Series)):
        acsm = pd.merge(acsm,vol, how='left', left_index=True,right_index=True)
        acsm['Vol PM1/PM2.5'] =acsm['V_SMPS[nm3/cm3]']/acsm['V_tot [nm3/cm3]']
    acsm['rho_tot [g/m3]'] = acsm['PM2.5 [ug/m3]']/acsm['V_tot [nm3/cm3]']*(10**9) # total density
    if (isinstance(vol, pd.Series)):
        acsm['PM2.5 [ug/m3]'] = acsm['PM2.5 [ug/m3]']
        acsm['PM1 [ug/m3]'] = acsm['V_SMPS[nm3/cm3]']/acsm['rho_tot [g/m3]']
        acsm['Mass PM1/PM2.5'] =acsm['PM1 [ug/m3]']/acsm['PM2.5 [ug/m3]']
    acsm['Nuetralization'] = (acsm['NH4[ug/m3]']/MW['NH4'])/(2*acsm['SO4[ug/m3]']/MW['SO4'] + acsm['NO3[ug/m3]']/MW['NO3'] + acsm['Chl[ug/m3]']/MW['Chl'])
    sulfate_kappa = []
    nitrate_kappa = []
    for N in acsm['Nuetralization'].to_numpy():
        if N<0.75:
            sulfate_kappa.append(0.8) # sulfuric acid
            nitrate_kappa.append(1) # more acidic nitrates
        elif N>1.25:
            sulfate_kappa.append(0.5) # ammonium bisulfate
            nitrate_kappa.append(0.67) # ammonium nitrate 
        else: 
            sulfate_kappa.append(0.6) # ammonium sulfate 
            nitrate_kappa.append(0.67) # ammonium nitrate
    acsm['k_SO4'] = sulfate_kappa
    acsm['k_NO3'] = nitrate_kappa
    acsm['k_total'] = (acsm['k_org']*acsm['V_org [nm3/cm3]'] + acsm['k_SO4']*acsm['V_SO4 [nm3/cm3]'] + acsm['k_NO3']*acsm['V_NO3 [nm3/cm3]'])/acsm['V_tot [nm3/cm3]']
    acsm[acsm<0] = 0
    # with pd.option_context('display.max_rows', 5, 'display.max_columns', None):
    #     input(acsm)
    return acsm


f2 = [r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_combined_PM25_2024.csv",
      r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_combined_PM25_2025.csv",
      r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_combined_PM25_2026.csv"]
fbc = [r"C:\Users\bensy\Documents\Research\AE33_Data\Aethelometer_Corrected_1hr_STP_2024.csv",
       r"C:\Users\bensy\Documents\Research\AE33_Data\Aethelometer_Corrected_1hr_STP_2025.csv",
       r"C:\Users\bensy\Documents\Research\AE33_Data\Aethelometer_Corrected_1hr_STP_2026.csv"]
vol =[r"C:\Users\bensy\Documents\Research\SMPS_2026\2026_SMPS_VolumeSizeDist_1hr_clean_stp.csv",
      r"C:\Users\bensy\Documents\Research\SMPS_2024\2024_SMPS_VolumeSizeDist_1hr_clean_stp.csv",
      r"C:\Users\bensy\Documents\Research\SMPS_2025\2025_SMPS_VolumeSizeDist_1hr_clean_stp.csv"]
vol_d= smps_vol(vol)
vol_d.rename('V_SMPS[nm3/cm3]',inplace=True)
print(vol_d.head)
pm25,specs = PM25_data(f2)
ae33 = AE33_data(fbc)
BC = ae33['BC_880nm_corr(ug/m3)']
BC.rename('BC[ug/m3]',inplace=True)
BC.dropna(how='all', inplace=True)
# input(BC)
# monthly_box_call(pm25,'PM2.5 [ug/m3]',y_label=r'Mass[ug/cm3]',title = 'Monthly AQS Aerosol Mass')
f = r"C:\Users\bensy\Documents\Research\20260512_ACSM_240529-260511_speciated_allmz(1).xlsx"
acsm = ACSM_data(f,BC = BC,vol=vol_d,freq='h')
acsm = acsm.dropna(axis='rows',how='all')
acsm.to_csv(r"C:\Users\bensy\Documents\Research\ACSMdata_240529-260511_incBC.csv")
vol_d.rename('V_tot [nm3/cm3]',inplace=True)
print(vol_d.to_frame()['V_tot [nm3/cm3]'])
print(acsm['V_tot [nm3/cm3]'])
mass_ratio = acsm['Mass PM1/PM2.5'].to_frame().dropna()
monthly_box_call(mass_ratio,'Mass PM1/PM2.5',y_label=r'pm1/pm2.5',title = 'Monthly Mass ratio between PM1 and PM2.5')
# monthly_box_call(acsm,'k_org',y_label=r'k$_{pm2.5}$',title = 'Monthly Organic Hygroscopicity Data')
monthly_box_call([vol_d.to_frame(),acsm],'V_tot [nm3/cm3]',y_label=r'V$_{tot}$ [nm3/cm3]',title = 'Monthly Aerosol volume',keys=['SMPS','ACSM'])