"""
Date: 6/10/2026
Author: Ben Sykes
Purpose: Process through raw CCN data and convert it according to ACTRIS formating standard
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from os.path import expanduser 
from CCN_process import *
try:
    from CCN_EBAS_convert import ebas_genfile
    ebas = True
except: 
    print('CCN_EBAS_convert script not accesable -> EBAS file not generated')
    ebas = False

'''EXAMPLE SET UP FOR EASE OF USE'''
'''Your exact folder structure will almost certainly look
different than this, but this is the general values you
will need to provide to run this processing in 'dev' 
mode. This process is also self constrained to run with
user inputs alone rather than supplied inputs in the 
global environment.'''
# turning 'dev' mode on. Turns off most user inputs, but assumes constants are provided.
dev = True
# Set up files in. My yearly cleaned files are all in the same folder with the structure "CCN_Clean_#YEAR#_1min.csv"
year = '2026'
file = expanduser(f"~/Documents/Research/CCN/CCN_Clean_{year}_1min.csv")                    # <- Necessary Input 
if year =='2026':
    ini_file = expanduser("~/Documents/Research/CCN 100 26.ini")                            # <- Necessary Input 
else:
    ini_file = expanduser("~/Documents/Research/CCN 100.ini") 
# set up files/folders out for levels 0-2 and the EBAS formatted data
file_out_lvl2 = expanduser(f"~/Documents/Research/CCN/Processed/CCN_lvl2_{year}_1hr.csv")   # <- Necessary Input 
file_out_lvl1 = expanduser(f"~/Documents/Research/CCN/Processed/CCN_lvl1_{year}_1hr.csv")   # <- Necessary Input 
file_out_lvl0 = expanduser(f"~/Documents/Research/CCN/Processed/CCN_lvl0_{year}_1hr.csv")   # <- Necessary Input 
ebas_out = expanduser("~/Documents/Research/CCN/Processed/")                                # <- Necessary Input if EBAS functionality provided
# provide known bad datas to remove from processing 
bad_dates = [[pd.to_datetime('10/01/2025 00:00:00'),pd.to_datetime('12/01/2025 00:00:00')], # <- Optional Input
             [pd.to_datetime('01/01/2026 00:00:00'),pd.to_datetime('01/09/2026 18:00:00')]]


def main():
    #region 
    '''====1 Read In Files 1====+++
    Read in the minute resolved files and CCN.ini files. If dev mode
    is active it assumes files paths were provided within the global 
    environment in an attempt to minimize required inputs.    
    +++====1 Read In Files 1===='''
    global file, ini_file, file_out, ebas_out, bad_dates, ebas
    if not dev:
        file = input('Provide path for CCN yearly csv file...')
        ini_file = input('Provide path for CCN.ini csv file...')
    df, dt_dct= readin(file)                                        # read in minute resolved data to DataFrame
    TGdum, slope_i, intercept_i, ss_list, date = readini(ini_file)  # read in init file to Metadata DataFrame
    ss_vals = [0.1,0.15,0.25,0.4,0.7]                               # set expected ss values 
    #endregion 
    #region 
    '''====2 Calculate SS and Apply Minute Flags 2====+++
    Calculate super saturation and Temperature Gradient using T1 and 
    T2 using the Khoeler curve assumption. 
    -- Minute Flags --
    I: ss_flag: denotes when the calculated ss% deviates by more than 20% off of the
                ss% machine set point
    II: Q_flag: denotes when the sample flow deviates more than 5% off of the mean
    III: T1_flag1: denotes when T1 is greater than 30*C
    IV: T1_flag2: denotes when T1 is less than T_inlet
    V: T1_flag3: denotes when T1 deviates more than 5*C from the mean T1 value
    +++====2 Calculate SS 2===='''
    df['ss_slope'],df['ss_intercept'] = slope_i, intercept_i
    T1= df['T1(C)'].to_numpy()
    T2= df['T2(C)'].to_numpy()
    slope = df['ss_slope'].to_numpy()
    intercept = df['ss_intercept'].to_numpy()
    # calculate super saturation and temperature gradient
    df['TG(C)_calc'],df['ss(%)_calc'] = sup_sat(T1,T2,A=slope, B= intercept)
    # ss flag calculation for the weighted correction
    # Flag bad datasets
    df['ss_dev'] = 100*np.abs(df['ss(%)_calc'].to_numpy() - df['ss(%)_setpt'].to_numpy())/((df['ss(%)_calc'].to_numpy() + df['ss(%)_setpt'].to_numpy())/2)
    df['ss_flag'] = (df['ss_dev'].to_numpy()>20.0).astype(int)
    Q1= df['Q(lpm)_sample'].to_numpy()
    df['Q_flag'] = (np.abs((Q1-np.nanmean(Q1))/np.nanmean(Q1)*100)>5).astype(int)                           #sample flow should stay within about 5%
    df['T1_flag1'] = (df['T1(C)'].to_numpy()>30.0).astype(int)                                              #T1 less than 30 deg
    df['T1_flag2'] = (df['T1(C)'].to_numpy()>df['T(C)_inlet'].to_numpy()).astype(int)                       #T1 less than Tinlet 
    df['T1_flag3'] = (np.abs(df['T1(C)'].to_numpy() - np.nanmean(df['T1(C)'].to_numpy()))>5).astype(int)    #T1 should stay within a 10ish degree band
    #endregion 
    #region
    '''====3 Apply Corrections 3====+++
    I: Calculate hourly averages by averaging together first by ss% set points
        and then by averaging over the hour
    II: Calculate the weighted correction by linear correcting the concentration
        values from the ss% set point, to the calculated ss% value
    III: Calculate the STP correction by adjusting flow values from ATP to STP
        and keep the STP set points in the Dataframe
    +++====3 Apply Corrections 3===='''
    #Keep the only time averaged df and the fully corrected df 
    ccn_corr_cols = ['N(cm-3)_cor_setpt0.1','N(cm-3)_cor_setpt0.15','N(cm-3)_cor_setpt0.25','N(cm-3)_cor_setpt0.4','N(cm-3)_cor_setpt0.7']
    df0 = time_avg_ss(df.copy(),ss_vals=ss_vals)                        # This version is only time-averaged
    df2 = time_avg_ss(df.copy(),ss_vals=ss_vals,clean = True)           # This version is fully cleaned at the minute resolution
    df1= weighted_corr(df0.copy(), ss_vals=ss_vals, param='N(cm-3)')    # Apply the ss correction to concentration
    df1 = stp_corr(df1, ccn_corr_cols,)                                 # Correct to STP
    df2 = weighted_corr(df2.copy(), ss_vals=ss_vals, param='N(cm-3)')   # Apply the ss correction to concentration
    df2 = stp_corr(df2, ccn_corr_cols,)                                 # Correct to STP
    #Drop unnecesary columns from all outputs 
    df0.dropna(axis =1, how='all', inplace=True)
    df0.dropna(axis =0,thresh = 5, inplace = True)
    df0 = df0.drop(columns=['T2(C)', 'T3(C)','TG(C)_calc','N(cm-3)','ss(%)_setpt', 'TG(C)_setpt','ss(%)_calc', 'ss_dev']) 

    df1.dropna(axis =1, how='all', inplace=True)
    df1.dropna(axis =0,thresh = 5, inplace = True)
    df1= df1.drop(columns=['T2(C)', 'T3(C)','TG(C)_calc','N(cm-3)','ss(%)_setpt', 'TG(C)_setpt','ss(%)_calc', 'ss_dev']) 

    df2.dropna(axis =1, how='all', inplace=True)
    df2.dropna(axis =0,thresh = 5, inplace = True)
    df2 = df2.drop(columns=['T2(C)', 'T3(C)','TG(C)_calc','N(cm-3)','ss(%)_setpt', 'TG(C)_setpt','ss(%)_calc', 'ss_dev']) 
    #endregion 
    #region 
    '''====4 Apply Hourly Flags 4====+++
    Apply the following QA Hourly flags for the CCN Data:
    -- Hourly Flags --
    I: integrity_flag: denotes when more than 25% of the hour (15 mins) is missing 
        from the averages
    II: N_flag: denotes when concentration of any ss% setpoint is greater than 5000#/cm^3
    ------
    All flags are LOW when the data is behaving normally and go HIGH when their conditions 
    are met. Flags are combined into a boolean flag code with flag I being the MSB and
    Flag VII being the LSB which is converted to a decimal equivalent for the csv output. 
    Minute resolved flags are used as an hourly flag in df0 and df1 if >50% of the hour 
    contained a flag. 
    +++====4 Apply Flags 4===='''
    # Level 0 
    for c in df0.columns.to_numpy():
        if ('ss(%)_calc_setpt' in c)|('(C)' in c)|('(hPA)' in c)|('(lpm)' in c):
            val = df0.pop(c)
            df0[c] = val
    ss_slope1 = df0.pop('ss_slope')
    ss_int1 = df0.pop('ss_intercept')
    df0['integrity_flag'] = (df0['avg_complete'].to_numpy()<.75).astype(int)
    df0['N_flag'] = ((df0['N(cm-3)_avg_setpt0.7'].to_numpy()>5000.0)|(df0['N(cm-3)_avg_setpt0.4'].to_numpy()>5000.0)|(df0['N(cm-3)_avg_setpt0.25'].to_numpy()>5000.0)|(df0['N(cm-3)_avg_setpt0.15'].to_numpy()>5000.0)|(df0['N(cm-3)_avg_setpt0.1'].to_numpy()>5000.0)).astype(int) #concentration less than 5000
    for c in ['ss_flag','integrity_flag','Q_flag','N_flag','T1_flag1','T1_flag2','T1_flag3']:
        val = df0.pop(c)
        df0[c] = np.round(val)
    flags = df0[['ss_flag','integrity_flag','Q_flag','N_flag','T1_flag1','T1_flag2','T1_flag3']].to_numpy()
    flag_code = [int(s, base=2) for s in np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags)] # generate a binary representation of the flags
    df0['flag_code'] = flag_code
    df0['date_run'] = pd.to_datetime('now',utc=True).date()
    df0['date_param'] = pd.to_datetime(date, utc=True).date()
    df0['ss_slope'] =ss_slope1
    df0['ss_int'] = ss_int1
    try:
        df1.pop('check')
    except:
        'do nothing'
    # Drop rows where all concentrations are 0 or NaN
    conc = [col for col in df0.columns.to_numpy() if 'N(cm-3)' in col]
    mask = ((df0[conc] == 0) | (df0[conc].isna())).all(axis=1)
    df0 = df0.loc[~mask]

    # Level 1 
    for c in df1.columns.to_numpy():
        if ('ss(%)_calc_setpt' in c)|('(C)' in c)|('(hPa)' in c)|('(lpm)' in c):
            val = df1.pop(c)
            df1[c] = val
    T_act = df1.pop('T(C)_sample')
    T_stp = df1.pop('T(C)_stp')
    P_act = df1.pop('P(hPA)_sample')
    P_stp = df1.pop('P(hPA)_stp')
    df1['T(C)_stp'] = T_stp
    df1['T(C)_sample'] = T_act
    df1['P(hPa)_stp'] = P_stp
    df1['P(hPa)_sample'] = P_act
    avg_comp = df1.pop('avg_complete')
    ss_slope = df1.pop('ss_slope')
    ss_int = df1.pop('ss_intercept')
    df1['avg_complete'] = avg_comp
    df1['integrity_flag'] = (df1['avg_complete'].to_numpy()<.75).astype(int)
    df1['N_flag'] = ((df1['N(cm-3)_cor_stp_setpt0.7'].to_numpy()>5000.0)|(df1['N(cm-3)_cor_stp_setpt0.4'].to_numpy()>5000.0)|(df1['N(cm-3)_cor_stp_setpt0.25'].to_numpy()>5000.0)|(df1['N(cm-3)_cor_stp_setpt0.15'].to_numpy()>5000.0)|(df1['N(cm-3)_cor_stp_setpt0.1'].to_numpy()>5000.0)).astype(int) #concentration less than 5000
    for c in ['ss_flag','integrity_flag','Q_flag','N_flag','T1_flag1','T1_flag2','T1_flag3']:
        val = df1.pop(c)
        df1[c] = np.round(val)
    flags = df1[['ss_flag','integrity_flag','Q_flag','N_flag','T1_flag1','T1_flag2','T1_flag3']].to_numpy()
    flag_code = [int(s, base=2) for s in np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags)]
    df1['flag_code'] = flag_code
    df1['date_run'] = pd.to_datetime('now',utc=True).date()
    df1['date_param'] = pd.to_datetime(date, utc=True).date()
    df1['ss_slope'] = ss_slope
    df1['ss_int'] = ss_int
    try:
        df1.pop('check')
    except:
        'do nothing'
    # Drop rows where all concentrations are 0 or NaN
    conc = [col for col in df1.columns.to_numpy() if 'N(cm-3)' in col]
    mask = ((df1[conc] == 0) | (df1[conc].isna())).all(axis=1)
    df1 = df1.loc[~mask]
        
    # Level 2
    for c in df2.columns.to_numpy():
        if ('ss(%)_calc_setpt' in c)|('(C)' in c)|('(hPa)' in c)|('(lpm)' in c):
            val = df2.pop(c)
            df2[c] = val
    T_act = df2.pop('T(C)_sample')
    T_stp = df2.pop('T(C)_stp')
    P_act = df2.pop('P(hPA)_sample')
    P_stp = df2.pop('P(hPA)_stp')
    df2['T(C)_stp'] = T_stp
    df2['T(C)_sample'] = T_act
    df2['P(hPa)_stp'] = P_stp
    df2['P(hPa)_sample'] = P_act
    avg_comp = df2.pop('avg_complete')
    ss_slope = df2.pop('ss_slope')
    ss_int = df2.pop('ss_intercept')
    df2['avg_complete'] = avg_comp
    df2['integrity_flag'] = (df2['avg_complete'].to_numpy()<.75).astype(int)
    df2['N_flag'] = ((df2['N(cm-3)_cor_stp_setpt0.7'].to_numpy()>5000.0)|(df2['N(cm-3)_cor_stp_setpt0.4'].to_numpy()>5000.0)|(df2['N(cm-3)_cor_stp_setpt0.25'].to_numpy()>5000.0)|(df2['N(cm-3)_cor_stp_setpt0.15'].to_numpy()>5000.0)|(df2['N(cm-3)_cor_stp_setpt0.1'].to_numpy()>5000.0)).astype(int) #concentration less than 5000
    for c in ['ss_flag','Q_flag','T1_flag1','T1_flag2','T1_flag3']:
        val = df2.pop(c) # drop minute collected flags since no longer needed
    flags = df2[['integrity_flag','N_flag']].to_numpy()
    flag_code = [int(s, base=2) for s in np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags)]
    df2['flag_code'] = flag_code
    df2['date_run'] = pd.to_datetime('now',utc=True).date()
    df2['date_param'] = pd.to_datetime(date, utc=True).date()
    df2['ss_slope'] = ss_slope
    df2['ss_int'] = ss_int
    try:
        df2.pop('check')
    except:
        'do nothing'
    # Drop rows where all concentrations are 0 or NaN
    conc = [col for col in df2.columns.to_numpy() if 'N(cm-3)' in col]
    mask = ((df2[conc] == 0) | (df2[conc].isna())).all(axis=1)
    for date in bad_dates: # remove periods of suspected bad data
        mask |= (df2.index >= date[0]) & (df2.index <= date[-1])
    df2 = df2.loc[~mask]

    #endregion 
    '''====5 Generate Outputs 5====+++
    Output the data to a ready to use CSV and a NASA AMES files as desired by Actris.
    +++====5 Generate Outputs 5===='''
    #region 
    if not dev:
        file_out = input('Provide output path for processed CCN yearly file...')
        if ebas:
            ebas_out = input('Provide output folder for formated CCN EBAS file...')
    df2.to_csv(file_out_lvl2)
    df1.to_csv(file_out_lvl1)
    df0.to_csv(file_out_lvl0)
    if ebas:
        CCN_EBAS(file_out_lvl2, ebas_out, ss_vals)
    print(f"Finished, files generated at {ebas_out}")
    #endregion 

if __name__:
    main()
