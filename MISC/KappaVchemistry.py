"""
Date: 4/25/2026
Author: Ben Sykes
Purpose: Compare kappa with 
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import least_squares as LSfit
from datetime import datetime 
# from ACSM import ACSM_data
from plotgen import monthly_box_call, hourly_box_call,var_box_call,line_plot, scat_call
from scipy.optimize import brentq
import time
import sys
import re
from scipy.interpolate import interp1d
from scipy.special import erf
pd.set_option('mode.chained_assignment', None)
plt.rcParams['font.size'] = 15
plt.rcParams['axes.titlesize'] =18

if __name__ == '__main__':
    # chem_file = r"C:\Users\bensy\Documents\Research\ACSMdata_240529-260511.csv"
    chemBC_file = r"C:\Users\bensy\Documents\Research\ACSMdata_240529-260511_incBC.csv"
    # acsm = ACSM_data(chem_file, freq='h')
    # acsm = pd.read_csv(chem_file, index_col=0)
    # acsm.index = pd.to_datetime(acsm.index, format='mixed')
    acsm = pd.read_csv(chemBC_file, index_col=0)
    acsm.index = pd.to_datetime(acsm.index, format='mixed')
    kappa_file = r"C:\Users\bensy\Documents\Research\Dcrit_Kappa_calculations.csv"
    kappa = pd.read_csv(kappa_file, index_col=0)
    kappa.index = pd.to_datetime(kappa.index, format='mixed')
    for col in kappa.columns.to_numpy():
        if 'N(cm-3)' in col:
            ss = list(col.split('setpt'))[-1]
            kappa[f'Fact(ss={ss})'] = kappa[col]/kappa['Total Concentration (#/cm³)']
    Data_check = pd.merge(acsm,kappa,left_index = True, right_index = True, how='right')
    Data_check.index.rename('Datetime', inplace=True)

    Data_check.to_csv(r"C:\Users\bensy\Documents\Research\CCN_SMPS_chem.csv")
    kVcols= [col for col in acsm.columns.to_numpy() if ('V_' in col)|('k' in col)|('O/C'in col)|('H/C'in col)|('HOA' in col)|('OOA' in col)]
    kV= acsm.copy()[kVcols]
    kV = kV.dropna(how='all',axis= 'rows')
    # kVcols= [col for col in acsm.columns.to_numpy() if ('V_' in col)|('k' in col)]
    # kBC= acsm.copy()[kVcols]
    # kBC = kBC.dropna(how='all',axis= 'rows')
    # kBC = kBC.add_suffix('BC')
    # # kBCcols = [col for col in kBC.columns.to_numpy() if 'k' in col]
    # kBC = kBC[kBCcols]
    # acsm=acsm.resample('MS').median()

    # input(kV)
    freq = 'd'
    print(Data_check)
    stds = kappa.copy().resample(freq).std()
    stds = stds.add_suffix(' std')
    means = kappa.copy().resample(freq).median()
    kcols = [col for col in means.columns.to_numpy() if ('Kappa' in col)|('Fact' in col)|('N' in col)]
    k = kappa.copy()[kcols]
    freq_c = 'h'
    kV = pd.merge(kV.copy().resample(freq_c).median(),k.copy().resample(freq_c).median(),left_index = True, right_index = True)
    kACSM = kV[kVcols]
    kCCN = kV[kcols]
    kV = kV.dropna(how='any',axis= 'rows')
    kept_months = kV.index.month.to_numpy(dtype=str)
    kept_years = kV.index.year.to_numpy(dtype=str)
    kept_dates = ["/".join(i) for i in zip(kept_months,kept_years)]
    # input(kept_dates)
    # input(kV.columns)
    kV.columns = kV.columns.str.replace('Kappa', 'k')
    kV['MAF'] = 1 - kV['V_BC [nm3/cm3]']/kV['V_tot [nm3/cm3]']
    # input(kV.columns)
    kV['Δk(ACSM-ss0.15)']= kV['k_total']-kV['k(ss=0.15)']
    # kV['Δk(ACSM&BC-ss0.15)']= kV['k_totalBC']-kV['k(ss=0.15)']
    # kV['Δk(ACSM-ACSM&BC)']= kV['k_total']-kV['k_totalBC']
    kCCN.columns = kCCN.columns.str.replace('Kappa', 'k')
    kV.columns = kV.columns.str.replace(' Est', '')
    kCCN['k_total'] = kV['k(ss=0.15)']
    kV['k_org(ss=0.15)'] = ((kV['k(ss=0.15)']*kV['V_tot [nm3/cm3]'] -kV['k_SO4']*kV['V_SO4 [nm3/cm3]']-kV['k_NO3']*kV['V_NO3 [nm3/cm3]'])/
                            (kV['V_tot [nm3/cm3]'] -kV['V_SO4 [nm3/cm3]']-kV['V_NO3 [nm3/cm3]']))
    kV['k_inorg(ss=0.25)'] = (kV['k(ss=0.25)']*kV['V_tot [nm3/cm3]'] -kV['k_org']*kV['V_org [nm3/cm3]'])/(kV['V_tot [nm3/cm3]'] -kV['V_org [nm3/cm3]'])
    kV['k_inorg(ss=0.15)'] = (kV['k(ss=0.15)']*kV['V_tot [nm3/cm3]'] -kV['k_org']*kV['V_org [nm3/cm3]'])/(kV['V_tot [nm3/cm3]'] -kV['V_org [nm3/cm3]'])

    mask = (kV.index <= pd.to_datetime('1/01/2025'))#werid data from 2024
    kV = kV[~mask]
    kV.to_csv(r"C:\Users\bensy\Documents\Research\kappa_Chem.csv")
    ee_file = r"C:\Users\bensy\Documents\Research\EarthEngine\MODIS_monthly_NDVI_FPAR_2024-06-01-2026-07-14_350.0km.csv"
    smoke_file = r"C:\Users\bensy\Documents\Research\SmokeData\dailySmoke.csv"
    smoke = pd.read_csv(smoke_file, index_col=0)
    smoke.index = pd.to_datetime(smoke.index, yearfirst=True)
    # input(smoke)
    kV = pd.merge(kV.copy().resample(freq_c).median(),smoke.copy(),left_index = True, right_index = True)
    ee = pd.read_csv(ee_file, index_col=0)
    ee.index = pd.to_datetime(ee.index, yearfirst=True)
    ee = ee.copy().resample('MS').median()
    print(ee)
    Data_out = pd.merge(acsm,means,left_index = True, right_index = True,how ='right')
    # Data_out.to_csv(r"C:\Users\bensy\Documents\Research\CCNEnviroCond.csv")
    Data_out = pd.merge(stds,Data_out,left_index = True, right_index = True)
    Data_out = Data_out.resample('MS').mean()
    print(Data_out)
    Data_out = pd.merge(Data_out, ee, left_index = True, right_index = True)
    kV.to_csv(r"C:\Users\bensy\Documents\Research\CCNEnviroCond.csv")
    cdf = Data_out.copy()
    print(cdf)
    cdf['Year'] = cdf.index.year.to_numpy()
    cdf['Month'] = cdf.index.month.to_numpy()
    normed =(cdf-cdf.min())/(cdf.max()-cdf.min())
    normed = normed.dropna(how='all')
    cdf = cdf[(cdf['Year'] == 2025)|(cdf['Year'] == 2026)]
    t5 = cdf[cdf['Year'] == 2025]
    t6 = cdf[cdf['Year'] == 2026]
    t5 = t5[t5['Month']<=t6['Month'].to_numpy()[-1]]
    # print(t5)
    # print(t6)
    avgs = cdf.groupby('Month').mean().reset_index()
    avgs = avgs[avgs['Month']<=t6['Month'].to_numpy()[-1]]
    print(avgs)
    for col in avgs: 
        avgs[f'yearly Δ{col}[%]'] = (t6[col].to_numpy() - t5[col].to_numpy())/np.nanmean(np.asarray([t5[col].to_numpy(),t6[col].to_numpy()])) *100
    avgs = avgs.set_index('Month')
    print(normed.index.to_numpy())
    print(normed.columns.to_numpy())
    print(kV.columns)
    # kV= kV.copy().resample('d').median()
    bad_date = [pd.to_datetime('10/01/2025'),pd.to_datetime('12/01/2025')]
    mask = (kV.index >= bad_date[0]) & (kV.index <= bad_date[-1])
    kV = kV[~mask]
    spring = [[pd.to_datetime('03/01/2025'),pd.to_datetime('06/01/2025')],[pd.to_datetime('03/01/2026'),pd.to_datetime('06/01/2026')]]
    summer = [[pd.to_datetime('06/01/2025'),pd.to_datetime('09/01/2025')],[pd.to_datetime('06/01/2026'),pd.to_datetime('09/01/2026')]]
    winter = [[pd.to_datetime('12/01/2024'),pd.to_datetime('03/01/2025')],[pd.to_datetime('12/01/2025'),pd.to_datetime('03/01/2026')]]
    sprmask = (kV.index <= pd.to_datetime('01/01/2024'))
    summask = (kV.index <= pd.to_datetime('01/01/2024'))
    winmask = (kV.index <= pd.to_datetime('01/01/2024'))
    spr26mask = (kV.index >= pd.to_datetime('03/01/2026')) & (kV.index <= pd.to_datetime('06/01/2026'))
    spr25mask = (kV.index >= pd.to_datetime('03/01/2025')) & (kV.index <= pd.to_datetime('06/01/2025'))
    for s in spring:
        sprmask |= (kV.index >= s[0]) & (kV.index <= s[-1])
    for s in summer:
        summask |= (kV.index >= s[0]) & (kV.index <= s[-1])
    for w in winter:
        winmask |= (kV.index >= w[0]) & (kV.index <= w[-1])

    # for ss in ['0.1','0.15','0.25','0.4','0.7']:
    #     kV[f'Fact(ss={ss})'] = acsm[f'Fact(ss={ss})']
    kV['ε_SO4'] = kV['V_SO4 [nm3/cm3]']/kV['V_tot [nm3/cm3]']
    kV['ε_org'] = kV['V_org [nm3/cm3]']/kV['V_tot [nm3/cm3]']
    kV['ε_NO3'] = kV['V_NO3 [nm3/cm3]']/kV['V_tot [nm3/cm3]']
    kV['ε_BC'] = kV['V_BC [nm3/cm3]']/kV['V_tot [nm3/cm3]']

    kV['Org Cont.'] = kV['k_org']/kV['k_total']
    kV['SO4 Contribution'] = kV['ε_SO4']*kV['k_SO4']
    kV['Fkorg'] = kV['k_org']/kV['k_total']
    kVSpring = kV.copy()[sprmask]
    kVsp26 = kV.copy()[spr26mask]
    kVsp25 = kV.copy()[spr25mask]
    kVSummer = kV.copy()[summask]
    kVWinter = kV.copy()[winmask]
    # input(kV.columns)
    kVsmoke = kV.copy()[kV['Density'] == 'Heavy']
    kVclear = kV.copy()[kV['Density'] == 'Light']
    var_box_call(kVSummer,'k(ss=0.7)','Density',y_label=r'k$_{CCN}$(ss=0.7)',title='Summer K Under Smokey Conditions')
    monthly_box_call([kVsmoke,kVclear],'k(ss=0.1)',y_label=r'ε$_{BC}$',title = 'Activation Under Smokey Conditions',keys=['Smoke Impacted','Clear'])
    kVsmoke = kVsmoke.add_suffix('sm',axis='columns')
    kVclear = kVclear.add_suffix('cl',axis='columns')
    kVSpring = kVSpring.add_suffix('Sp',axis='columns')
    kVSummer = kVSummer.add_suffix('Sm',axis='columns')
    kVWinter = kVWinter.add_suffix('Wt',axis='columns')

    # input(kV.columns)
    d_smoke = {'Fact(ss=0.1)sm':'ε_orgsm','Fact(ss=0.1)cl':'ε_orgcl',
                'Fact(ss=0.7)sm':'ε_orgsm','Fact(ss=0.7)cl':'ε_orgcl'}
    leg_smoke= {'Fact(ss=0.1)sm':r'Smoke(ss=0.1)','Fact(ss=0.1)cl':r'Clear(ss=0.1)','Fact(ss=0.7)sm':r'Smoke(ss=0.7)','Fact(ss=0.7)cl':r'Clear(ss=0.7)'}
    M,B,R = scat_call([kVsmoke,kVclear,kVsmoke,kVclear],d_smoke,x_label=r'ε$_{org}$',y_label=r'F$_{act}$', title=r'F$_{act}$ vs f$_{HOA}$', single=False,leg_dict=leg_smoke)#,xtrp=1.0)#
    
    d_season = {'Fact(ss=0.1)Wt':'HOA/OAWt','Fact(ss=0.1)Sm':'HOA/OASm',
                'Fact(ss=0.7)Wt':'HOA/OAWt','Fact(ss=0.7)Sm':'HOA/OASm'}
    # input(kVSummer)
    leg_season = {'Fact(ss=0.1)Wt':r'Winter(ss=0.1)','Fact(ss=0.1)Sm':r'Summer(ss=0.1)','Fact(ss=0.7)Wt':r'Winter(ss=0.7)','Fact(ss=0.7)Sm':r'Summer(ss=0.7)'}
    M,B,R = scat_call([kVWinter,kVSummer,kVWinter,kVSummer],d_season,x_label=r'HOA/OA$',y_label=r'F$_{act}$', title=r'F$_{act}$ vs f$_{HOA}$', single=False,leg_dict=leg_season)#,xtrp=1.0)#

    d_season = {'k(ss=0.1)':'ε_org','k_org':'ε_org'}
    leg_season = {'k(ss=0.1)':r'κ$_{CCN}$(ss=0.15)','k_org':r'κ$_{org}$'}
    M,B,R = scat_call([kV,kV],d_season,x_label=r'ε$_{org}$',y_label=r'κ', title=r'κ vs ε$_{org}$', single=False,leg_dict=leg_season,xtrp=1.0,verbose=False)#
    
    d_season = {'k(ss=0.1)Sm':'ε_orgSm','k_orgSm':'ε_orgSm'}
    leg_season = {'k(ss=0.1)Sm':r'κ$_{CCN}$(ss=0.1)','k_orgSm':r'κ$_{org}$'}
    M,B,R = scat_call([kVSummer,kVSummer],d_season,x_label=r'ε$_{org}$',y_label=r'κ', title=r'κ vs ε$_{org}[Summer 2025]$', single=False,leg_dict=leg_season,xtrp=1.0,verbose=False)#

    d_season = {'k(ss=0.1)Wt':'ε_orgWt','k_orgWt':'ε_orgWt'}
    leg_season = {'k(ss=0.1)Wt':r'κ$_{CCN}$(ss=0.1)','k_orgWt':r'κ$_{org}$'}
    M,B,R = scat_call([kVWinter,kVWinter],d_season,x_label=r'ε$_{org}$',y_label=r'κ', title=r'κ vs ε$_{org}[Winter]$', single=False,leg_dict=leg_season,xtrp=1.0,verbose=False)#

    # input(kVsp26.columns)
    d_season = {'k(ss=0.1)Sp':'ε_orgSp','k_orgSp':'ε_orgSp'}
    leg_season = {'k(ss=0.1)Sp':r'κ$_{CCN}$(ss=0.1)','k_orgSp':r'κ$_{org}$'}
    M,B,R = scat_call([kVSpring,kVSpring],d_season,x_label=r'ε$_{org}$',y_label=r'κ', title=r'κ vs ε$_{org}[Spring]$', single=False,leg_dict=leg_season,xtrp=1.0,verbose=False)#


    kVsp26 = kVsp26.add_suffix('26',axis='columns')
    kVsp25= kVsp25.add_suffix('25',axis='columns')
    d_season = {'k(ss=0.1)25':'ε_org25','k(ss=0.1)26':'ε_org26'}
    leg_season = {'k(ss=0.1)25':r'Spring 25','k(ss=0.1)26':r'Spring 26'}
    M,B,R = scat_call([kVsp25,kVsp26],d_season,x_label=r'ε$_{org}$',y_label=r'κ$_{CCN}$', title=r'κ$_{CCN}$(ss=0.1%) vs ε$_{org}$', single=False,leg_dict=leg_season,xtrp=1.0)#

    d_kappa = {'k(ss=0.1)':'ε_org','k(ss=0.4)':'ε_org','k(ss=0.7)':'ε_org'}#'ε_org''Fkorg'
    leg_dict = {'k(ss=0.1)':r'ss=0.1','k(ss=0.4)':r'ss=0.4','k(ss=0.7)':r'ss=0.7'}
    M,B,R = scat_call(kV,d_kappa,x_label=r'ε$_{BC}$',y_label='Kappa', title=r'κ$_{CCN}$ vs ε$_{BC}$', single=False,leg_dict=leg_dict) #xtrp=1.0,

    monthly_box_call(kV,'Fact(ss=0.1)',y_label=r'MAF estimate',title = 'Hygroscopic fraction')
    # monthly_box_call([kCCN,kACSM],'k_total',y_label=r'Hygroscopcity',title = 'Hygroscopicity from CCN and ACSM',keys=['CCN(ss=0.15)','ACSM'])
    monthly_box_call([kCCN,kACSM],'k_total',y_label=r'Hygroscopcity',title = 'Hygroscopicity from CCN and ACSM',keys=['CCN(ss=0.1)','ACSM+BC'])
    

    plt.ion()
    fig, ax = plt.subplots()
    ax.plot(kV['Δk(ACSM-ss0.15)'],'.-', label = r'ACSM-ss0.15')
    ax.legend(loc='center left',fontsize='medium')
    ax.set_ylabel('Δk')
    ax.set_xlabel('Date')
    ax.set_title('Deviation in Hygroscopcity between measurements')
    input('Press enter to exit plot...')
    ax.cla()

    fig, ax = plt.subplots()
    kV = kV.resample('ME').median()
    ax.plot(kV['k(ss=0.15)'], label ='CCN(ss=0.25)')
    ax.plot(kV['k_totalBC'], label = 'ACSM+BC')
    ax.legend()
    ax.set_ylabel('Hygroscopicity')
    ax.set_xlabel('Date')
    ax.set_title('Hygroscopicity Between Measurements')
    input('Press enter to exit plot...')
    ax.cla()

    fig, ax = plt.subplots()
    res = linregress(avgs['yearly ΔFact(ss=0.7)[%]'].to_numpy(),avgs['yearly ΔNDVI[%]'].to_numpy())
    m,b,r = res.slope,res.intercept,res.rvalue**2
    L, = ax.plot(avgs['yearly ΔFact(ss=0.7)[%]'].to_numpy(),avgs['yearly ΔNDVI[%]'].to_numpy(), label = f'{m}x+{b} | R2 = {r}', marker = '*',ls = '')
    fit, = ax.plot(avgs['yearly ΔFact(ss=0.7)[%]'].to_numpy(), m*avgs['yearly ΔFact(ss=0.7)[%]'].to_numpy()+b, color = L.get_color(), alpha = 0.75)
    ax.legend()
    ax.set_xlabel('ΔFact(ss=0.7)[%]')
    ax.set_ylabel('ΔNDVI[%]')
    ax.set_title("Correlation of Yearly Change('25->'26)")
    input('Press enter to exit plot...')
    ax.cla()

    fig, ax = plt.subplots()
    res = linregress(avgs['yearly ΔKappa(ss=0.25)[%]'].to_numpy(),avgs['yearly ΔNDVI[%]'].to_numpy())
    m,b,r = res.slope,res.intercept,res.rvalue**2
    L, = ax.plot(avgs['yearly ΔKappa(ss=0.25)[%]'].to_numpy(),avgs['yearly ΔNDVI[%]'].to_numpy(), label = f'{m}x+{b} | R2 = {r}', marker = '*',ls = '')
    fit, = ax.plot(avgs['yearly ΔKappa(ss=0.25)[%]'].to_numpy(), m*avgs['yearly ΔKappa(ss=0.25)[%]'].to_numpy()+b, color = L.get_color(), alpha = 0.75)
    ax.legend()
    ax.set_xlabel('ΔKappa(ss=0.25)[%]')
    ax.set_ylabel('ΔNDVI[%]')
    ax.set_title("Correlation of Yearly Change('25->'26)")
    input('Press enter to exit plot...')
    ax.cla()

    plt.ioff()
    
    print('finished')