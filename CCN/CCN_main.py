"""
Date: 1/29/2026
Author: Ben Sykes
Purpose: Process through raw CCN data and convert it according to ACTRIS formating standard
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import os
from os.path import expanduser 
import datetime as dt 
import scipy as sp
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
from CCN_EBAS_convert import ebas_genfile
from pathlib import Path
import time
from CCN_process import *

dev = True
file = expanduser("~\Documents\Research\CCN_Clean_2025_1min.csv")
ini_file = expanduser("~/Documents/Research/CCN 100.ini")
file_out = expanduser("~/Documents/Research/CCN_Processed_2025_1hr.csv")
ebas_out = expanduser("~/Documents/Research")


def main():
    '''====1 Read In Files 1====+++
    Read in the minute averaged files and CCN.ini files. If dev mode
    is active it assumes files paths were provided within the global 
    environment in an attempt to minimize required inputs.    
    +++====1 Read In Files 1===='''
    global file, ini_file, file_out, ebas_out
    #region 
    if not dev:
        file = input('Provide path for CCN yearly file...')
        ini_file = input('Provide path for CCN.ini file...')
    df, dt_dct= readin(file)
    TGdum, slope_i, intercept_i, ss_list, date = readini(ini_file)
    ss_vals = [0.1,0.15,0.25,0.4,0.7]
    #endregion 
    '''====2 Calculate SS 2====+++
    Calculate super saturation and Temperature Gradient using T1 and 
    T2 using the Khoeler curve assumption. 
    +++====2 Calculate SS 2===='''
    #region 
    df['ss_slope'],df['ss_intercept'] = slope_i, intercept_i
    T1= df['T1(C)'].to_numpy()
    T2= df['T2(C)'].to_numpy()
    slope = df['ss_slope'].to_numpy()
    intercept = df['ss_intercept'].to_numpy()
    df['TG(C)_calc'],df['ss(%)_calc'] = sup_sat(T1,T2,A=slope, B= intercept)
    #ss flag calculation for the weighted correction
    #Flag when the calculated ss is greater than 20% different than the machine set point
    df['ss_dev'] = 100*np.abs(df['ss(%)_calc'].to_numpy() - df['ss(%)_setpt'].to_numpy())/((df['ss(%)_calc'].to_numpy() + df['ss(%)_setpt'].to_numpy())/2)
    df['ss_flag'] = (df['ss_dev'].to_numpy()>20.0).astype(int)
    #endregion 
    '''====3 Apply Corrections 3====+++
    I: Calculate hourly averages by averaging together first by ss% set points
        and then by averaging over the hour
    II: Calculate the weighted correction by linear correcting the concentration
        values from the ss% set point, to the calculated ss% value
    III: Calculate the STP correction by adjusting flow values from ATP to STP
        and keep the STP set points in the Dataframe
    +++====3 Apply Corrections 3===='''
    #region 
    ccn_corr_cols = ['N(cm-3)_cor_setpt0.1','N(cm-3)_cor_setpt0.15','N(cm-3)_cor_setpt0.25','N(cm-3)_cor_setpt0.4','N(cm-3)_cor_setpt0.7']
    df = time_avg_ss(df,ss_vals=ss_vals)
    df = weighted_corr(df, ss_vals=ss_vals, param='N(cm-3)')
    df = stp_corr(df, ccn_corr_cols,)
    #Drop unnecesary columns
    df.dropna(axis =1, how='all', inplace=True)
    df.dropna(axis =0,thresh = 5, inplace = True)
    df = df.drop(columns=['T2(C)', 'T3(C)','TG(C)_calc','N(cm-3)','ss(%)_setpt', 'TG(C)_setpt','ss(%)_calc', 'ss_dev']) 
    #endregion 
    '''====4 Apply Flags 4====+++
    Apply the following QA flags for the CCN Data:
    I: ss_flag: denotes when the calculated ss% deviates by more than 20% off of the
        ss% machine set point
    II: integrity_flag: denotes when more than 25% of the hour (15 mins) is missing 
        from the averages
    III: Q_flag: denotes when the sample flow deviates more than 5% off of the mean
    IV: N_flag: denotes when concentration of any ss% setpoint is greater than 5000#/cm^3
    V: T1_flag1: denotes when T1 is greater than 30*C
    VI: T1_flag2: denotes when T1 is less than T_inlet
    VII: T1_flag3: denotes when T1 deviates more than 5*C from the mean T1 value
    ------
    All flags are LOW when the data is behaving normally and go HIGH when their conditions 
    are met. Flags are combined into a boolean flag code with flag I being the MSB and
    Flag VII being the LSB which is converted to a decimal equivalent for the csv output.
    +++====4 Apply Flags 4===='''
    #region 
    ss_flag = df.pop('ss_flag')# move all flags to end of dataframe
    Q = df['Q(lpm)_sample'].to_numpy()
    df['ss_variation'] = ss_flag
    df['ss_flag'] = np.round(ss_flag)
    df['integrity_flag'] = (df['avg_complete'].to_numpy()<.75).astype(int)
    df['Q_flag'] = (np.abs((Q-np.nanmean(Q))/np.nanmean(Q)*100)>5).astype(int) #sample flow should stay within about 5%
    df['N_flag'] = ((df['N(cm-3)_cor_stp_setpt0.7'].to_numpy()>5000.0)|(df['N(cm-3)_cor_stp_setpt0.4'].to_numpy()>5000.0)|(df['N(cm-3)_cor_stp_setpt0.25'].to_numpy()>5000.0)|(df['N(cm-3)_cor_stp_setpt0.15'].to_numpy()>5000.0)|(df['N(cm-3)_cor_stp_setpt0.1'].to_numpy()>5000.0)).astype(int) #concentration less than 5000
    df['T1_flag1'] = (df['T1(C)'].to_numpy()>30.0).astype(int) #T1 less than 30 deg
    df['T1_flag2'] = (df['T1(C)'].to_numpy()>df['T(C)_inlet'].to_numpy()).astype(int) #T1 less than Tinlet 
    df['T1_flag3'] = (np.abs(df['T1(C)'].to_numpy() - np.nanmean(df['T1(C)'].to_numpy()))>5).astype(int) #T1 should stay within a 10ish degree band
    flags = df[['ss_flag','integrity_flag','Q_flag','N_flag','T1_flag1','T1_flag2','T1_flag3']].to_numpy()
    flag = np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags)
    flag_code = [int(s, base=2) for s in np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags)]
    df['flag_code'] = flag_code
    df['date_run'] = pd.to_datetime('now',utc=True).date()
    df['date_param'] = pd.to_datetime(date, utc=True).date()
    df['ss_slope'] = df['ss_slope'].to_numpy()
    df['ss_int'] = df['ss_intercept'].to_numpy()
    #endregion 
    '''====5 Generate Outputs 5====+++
    Output the data to a ready to use CSV and a NASA AMES files as desired by Actris.
    +++====5 Generate Outputs 5===='''
    #region 
    if not dev:
        file_out = input('Provide output path for processed CCN yearly file...')
        ebas_out = input('Provide output folder for formated CCN EBAS file...')
    df.to_csv(file_out)
    CCN_EBAS(file_out,ebas_out, ss_vals)
    print(f"Finished, files generated at {file_out} and {ebas_out}")
    #endregion 

if __name__:
    main()



"""archived useful code"""
## slope for ss slope and intercept. Probably implemented in future iterations
# df, dt_dct= readin(file)
# TGdum, slope_i, intercept_i, ss_list, date = readini(ini_file)
# # slope_e, intercept_e = 17.06,0.73
# ss_vals = [0.1,0.15,0.25,0.4,0.7]
# start_date = pd.to_datetime('5/23/2025 00:00:00')
# end_date = pd.to_datetime('12/03/2025 00:00:00')
# timeD = pd.Timedelta(end_date-start_date).total_seconds()/60 #hours
# s_line = (slope_e -slope_i)/timeD*np.arange(0,timeD)+ slope_i
# i_line = (intercept_e -intercept_i)/timeD*np.arange(0,timeD)+ intercept_i
# ss_dict = {'Date String (YYYY-MM-DD hh:mm:ss) UTC':pd.date_range(start_date,end_date,freq='min',inclusive='left'), 'ss_slope':s_line, 'ss_intercept':i_line}
# ss_df = pd.DataFrame(ss_dict)
# ss_df =ss_df.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
# df =pd.merge(df, ss_df,left_index = True, right_index = True, how ='left')
# df[['ss_slope','ss_intercept']] = df[['ss_slope','ss_intercept']].bfill().ffill()
# input(df[['SS slope','SS intercept']])
