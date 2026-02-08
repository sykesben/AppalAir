"""
Date: 1/29/2026
Author: Ben Sykes
Purpose: Functions for processing CCN data
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

"""=====I Read in CCN File and CCN.ini file I====="""
def readin(path):
    """
    Takes in a path to a CCN file and generates outputs
    ----------
    Paramaters
    ++++++++++
    path : [str/path-like] Path to the CCN.ini file

    Returns
    ++++++++++
    data : [pandas DataFrame] read in data
    cols_rename : [dict] dictionary with verbose definition as the key and column name as the value
    """
    data = pd.read_csv(path,skiprows=lambda x:x==1)#read in csv skipping first row of verbose column headings
    data=data.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
    data.index.names = ['Datetime UTC']
    data.index = pd.to_datetime(data.index)
    cols_rename = {'particle number concentration (cm-3)': 'N(cm-3)','inlet temperature (°C)': 'T(C)_inlet','temperature of TEC 1 (°C)':'T1(C)','temperature of TEC 2 (°C)':'T2(C)'
                   ,'temperature of TEC 3 (°C)':'T3(C)','sample temperature (°C)': 'T(C)_sample','OPC temperature (°C)':'T(C)_OPC','nafion temperature (°C)':'T(C)_nafion',
                   'sample flow rate (lpm)': 'Q(lpm)_sample','sheath flow (lpm)':'Q(lpm)_sheath','reported supersaturation from onboard instrument calibration (%)':'ss(%)_setpt',
                   'sample pressure (hPa)':'P(hPA)_sample','temperature gradiant setpoint (°C)': 'TG(C)_setpt'}
    data.rename(columns=cols_rename,inplace=True)
    return data,cols_rename

def readini(path):
    """
    Takes in a path to a CCN.ini file and generates outputs
    ----------
    Paramaters
    ++++++++++
    path : [str/path-like] Path to the CCN.ini file

    Returns
    ++++++++++
    TGdum : [float] Temperature gradient offset for T3 calculation
    TGslope : [float] Intercept for Temperature Gradiant khoeler curve
    TGintercept : [float] Slope for Temperature Gradiant khoeler curve
    ss_list : [list] List of ss% set points
    date : [str] Date the ini file was generated
    """
    start = time.time()
    date = ''
    TGslope, TGintercept, TGdum = 0,0,0 #temperature gradient values
    date = ''
    ss_list = [0.1,0.15,0.25,0.4,0.7] #[0.1,0.15,0.2,0.25,0.4,0.6,0.7,0.8,1.0]
    file = open(path, "r")
    for line in file:
        #SS Table 1 0.SS Table Settings.SS %
        l = line.strip()
        if 'TG Dum' in l: 
            TGdum = float(list(l.split('='))[-1])
        elif 'Temp Gradient Slope' in l: 
            TGslope = float(list(l.split('='))[-1])
        elif 'Temp Gradient Y-intercept' in l: 
            TGintercept = float(list(l.split('='))[-1])
        elif 'Last Date Updated' in l: 
            date = (list(l.split('='))[-1])
    file.close()
    ss_list = np.unique(ss_list)
    time.sleep(1)
    end = time.time()
    print(f'Read in function time is {end-start}')
    return TGdum, TGslope, TGintercept, ss_list, date

"""=====II Calculate Super Saturation from Temp Gradient II====="""
def sup_sat(T1,T2,A=16.01,B=1.03):#A=16.01,B=1.03
    '''
    Super-Saturation Calculator from Temperature gradient. Assume linear relationship 
    between ss% and TG (for ss%>0.1%) such that TG = A*ss% + B.
    ----------
    Paramaters
    ++++++++++
    T1: [str] Temperature at TEC 1 [°C] (usually variable: 'T1_N12')
    T2: [str] Temperature at TEC 2 [°C] (usually variable: 'T2_N12')
    A: [float] Slope (default = 16.01)
    B: [float] Intercept (default = 1.03)

    Returns
    ++++++++++
    tg = temperature gradient len of T1 [C]
    ss: array of super saturation. Len of T1 [%]   
    '''
    tg=(T2-T1)*2 #Use this value according to actris formating, could use T3 instead of T2, but with TG adjustment
    ss = (tg - B)/A
    return tg,ss

"""=====III Average data over time and super saturation set point III====="""
def time_avg_ss(df, deltat='1h', ss_vals = [], ssflag = True):
    '''
    Groups values by machine set super saturation allowing for group time averaged values
    ----------
    Paramaters
    ++++++++++
    df : [pandas.Dataframe] Data to time average 
    deltat : [str] time period to re-average to (default: '1h', must be in '#n' format)
        possible values: m = minute, h = hour, d = day, M = month, Y = year
    ss_vals : [list or dict] possible values of ss from machine. If no value is given look at all recorded set points.
        IF a dictionary, the keys should be dates and the values ss lists to apply seperate ss values over
        different date ranges. (default: [])
    ss_flag: [list of bool] Does the measured ss devate more than 20% from the set point. (default: True)

    Returns
    +++++++
    df: [pandas.Dataframe] Time averaged data
    '''
    if isinstance(ss_vals, dict):
        start = time.time()
        dt_dict = {'m': 1, 'h': 60, 'D': 1440, 'M': 43830, 'Y': 525960}
        dt_num = float(''.join(filter(str.isdigit, deltat)))
        dt_minutes = dt_num * dt_dict[deltat[1]]
        df_new =pd.DataFrame()

        for date, ss_list in ss_vals.items():
            start, end = map(pd.to_datetime, date.split(' - '))
            data = df.loc[start:end]

            if len(ss_list) == 0: #if no machine super saturation set points, keep the unique supersaturations
                ss_list = np.unique(df['ss(%)_setpt'].to_numpy())
            
            # keep only stable ss setpoints
            data = data[data['ss(%)_setpt'] == data['ss(%)_setpt'].shift()]

            if ssflag:
                data = data[data['ss_flag'] == 0]

            pc_dict = {f'N(cm-3)_avg_setpt{s}': [] for s in ss_list}
            ss_dict = {f'ss(%)_calc_setpt{s}': [] for s in ss_list}
            tg_dict = {f'TG(C)_avg_setpt{s}': [] for s in ss_list}
            t1_dict = {f'T1(C)_avg_setpt{s}': [] for s in ss_list}
            t2_dict = {f'T2(C)_avg_setpt{s}': [] for s in ss_list}

            data_new = data.resample(deltat).mean() #Start at the top of the hour,day, month
            times = data_new.index.to_numpy()
            completeness = []

            for i in range(len(times)):
                ts = times[i]
                tf = ts +pd.Timedelta(dt_num, deltat[1])- pd.Timedelta(1,'min')
                slct = data.loc[ts:tf]
                completeness.append(len(slct)/dt) #how many minutes/vs minutes in an hour
                for ss in ss_list:
                    if ssflag:
                        pc = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)), 'N(cm-3)'].to_numpy())
                        cal = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'ss(%)_calc'].to_numpy())
                        tg = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'TG(C)_calc'].to_numpy())
                        t1 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'T1(C)'].to_numpy())
                        t2 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'T2(C)'].to_numpy())
                    else:
                        pc = np.nanmean(slct.loc[slct['ss(%)_setpt'] == ss, 'N(cm-3)'].to_numpy())
                        cal = np.nanmean(slct.loc[(slct['ss(%)_setpt'] == ss),'ss(%)_calc'].to_numpy())
                        tg = np.nanmean(slct.loc[slct['ss(%)_setpt'] == ss,'TG(C)_calc'].to_numpy())
                        t1 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)),'T1(C)'].to_numpy())
                        t2 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)),'T2(C)'].to_numpy())
                    pc_dict[f'N(cm-3)_avg_setpt{ss}'].append(pc)
                    ss_dict[f'ss(%)_calc_setpt{ss}'].append(cal)
                    tg_dict[f'TG(C)_avg_setpt{ss}'].append(tg)
                    t1_dict[f'T1(C)_avg_setpt{ss}'].append(t1)
                    t2_dict[f'T2(C)_avg_setpt{ss}'].append(t2)
            pc_cols = pd.DataFrame.from_dict(pc_dict)
            pc_cols.index = data_new.index.to_numpy()
            ss_cols = pd.DataFrame.from_dict(ss_dict)
            ss_cols.index = data_new.index.to_numpy()
            tg_cols = pd.DataFrame.from_dict(tg_dict)
            tg_cols.index = data_new.index.to_numpy()
            t1_cols = pd.DataFrame.from_dict(t1_dict)
            t1_cols.index = data_new.index.to_numpy()
            t2_cols = pd.DataFrame.from_dict(t2_dict)
            t2_cols.index = data_new.index.to_numpy()
            data_new = data_new.merge(pc_cols,right_index = True, left_index = True)
            data_new = data_new.merge(ss_cols,right_index = True, left_index = True)
            data_new = data_new.merge(tg_cols,right_index = True, left_index = True)
            data_new = data_new.merge(t1_cols,right_index = True, left_index = True)
            data_new = data_new.merge(t2_cols,right_index = True, left_index = True)
            data_new['avg_complete'] = completeness
            if df_new.empty:
                df_new = data_new
            else:
                df_new = pd.concat([df_new,data_new])
        end = time.time()
        print(f'time average function time is {end-start}')
        return df_new
    else: 
        start = time.time()
        ss_list= ss_vals
        dt_dict = {'m':1, 'h':60, 'D':1_440, 'M':43_830, 'Y': 525_960} # conversion factor to minutes
        dt_num = float(''.join([n for n in deltat if n.isdigit()]))
        dt = dt_num * dt_dict[deltat[1]]
        if len(ss_list) == 0:
            ss_list = np.unique(df['ss(%)_setpt'].to_numpy())
        pc_dict = {f'N(cm-3)_avg_setpt{s}': [] for s in ss_list}
        ss_dict = {f'ss(%)_calc_setpt{s}': [] for s in ss_list}
        tg_dict = {f'TG(C)_avg_setpt{s}': [] for s in ss_list}
        t1_dict = {f'T1(C)_avg_setpt{s}': [] for s in ss_list}
        t2_dict = {f'T2(C)_avg_setpt{s}': [] for s in ss_list}
        df_new = df.resample(deltat).mean() #Start at the top of the hour,day, month
        times = df_new.index.to_numpy()
        completeness = []
        for i in range(len(times)):
            ts = times[i]
            tf = ts +pd.Timedelta(dt_num, deltat[1])- pd.Timedelta(1,'min')
            slct = df.loc[ts:tf]
            completeness.append(len(slct)/dt) #how many minutes/vs minutes in an hour
            for ss in ss_list:
                if ssflag:
                    pc = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)), 'N(cm-3)'].to_numpy())
                    cal = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'ss(%)_calc'].to_numpy())
                    tg = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'TG(C)_calc'].to_numpy())
                    t1 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'T1(C)'].to_numpy())
                    t2 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0)),'T2(C)'].to_numpy())
                else:
                    pc = np.nanmean(slct.loc[slct['ss(%)_setpt'] == ss, 'N(cm-3)'].to_numpy())
                    cal = np.nanmean(slct.loc[(slct['ss(%)_setpt'] == ss),'ss(%)_calc'].to_numpy())
                    tg = np.nanmean(slct.loc[slct['ss(%)_setpt'] == ss,'TG(C)_calc'].to_numpy())
                    t1 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)),'T1(C)'].to_numpy())
                    t2 = np.nanmean(slct.loc[((slct['ss(%)_setpt'] == ss)),'T2(C)'].to_numpy())
                pc_dict[f'N(cm-3)_avg_setpt{ss}'].append(pc)
                ss_dict[f'ss(%)_calc_setpt{ss}'].append(cal)
                tg_dict[f'TG(C)_avg_setpt{ss}'].append(tg)
                t1_dict[f'T1(C)_avg_setpt{ss}'].append(t1)
                t2_dict[f'T2(C)_avg_setpt{ss}'].append(t2)
        pc_cols = pd.DataFrame.from_dict(pc_dict)
        pc_cols.index = df_new.index.to_numpy()
        ss_cols = pd.DataFrame.from_dict(ss_dict)
        ss_cols.index = df_new.index.to_numpy()
        tg_cols = pd.DataFrame.from_dict(tg_dict)
        tg_cols.index = df_new.index.to_numpy()
        t1_cols = pd.DataFrame.from_dict(t1_dict)
        t1_cols.index = df_new.index.to_numpy()
        t2_cols = pd.DataFrame.from_dict(t2_dict)
        t2_cols.index = df_new.index.to_numpy()
        df_new = df_new.merge(pc_cols,right_index = True, left_index = True)
        df_new = df_new.merge(ss_cols,right_index = True, left_index = True)
        df_new = df_new.merge(tg_cols,right_index = True, left_index = True)
        df_new = df_new.merge(t1_cols,right_index = True, left_index = True)
        df_new = df_new.merge(t2_cols,right_index = True, left_index = True)
        df_new['avg_complete'] = completeness
        end = time.time()
        print(f'Time average function time is {end-start}')
        return df_new

"""=====IV Apply a weighted linear correction to the CCN particle number to ss% IV====="""
def rowwise_linfit(X,X_new, Y):
    """
    Fits Y_i = a_i * X_i + b_i for each row i.
    ----------
    Paramaters
    ++++++++++
    X : [array-like] Independent variable values w/ shape (n_rows, n_points)
    Y : [array-like] Dependent variable values w/ shape (n_rows, n_points)
    
    Returns
    ++++++++++
    slopes : [ndarray] Fitted slopes for each row w/ shape (n_rows,)
    intercepts : [ndarray] Fitted intercepts for each row w/ shape (n_rows,)
    Y_fit : [ndarray] Fitted values using the row-wise linear models w/ shape (n_rows, n_points)
    """
    start = time.time()
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    X_new = np.asarray(X_new, dtype = float)

    # Compute means per row
    X_mean = X.mean(axis=1, keepdims=True)
    Y_mean = Y.mean(axis=1, keepdims=True)

    # Compute slope (a_i) and intercept (b_i)
    numerator = np.sum((X - X_mean) * (Y - Y_mean), axis=1)
    denominator = np.sum((X - X_mean)**2, axis=1)
    slopes = numerator / denominator
    intercepts = (Y_mean.squeeze() - slopes * X_mean.squeeze())

    # Compute fitted values
    Y_fit = slopes[:, None] * X_new + intercepts[:, None]
    Y_fit = np.maximum(Y_fit, 0)
    end = time.time()
    print(f'linear fit function time is {end-start}')
    return Y_fit

def weighted_corr(df,ss_vals,param):
    """
    Applies a weighted average to generate a linear fit for N vs ss
    ----------
    Paramaters
    ++++++++++
    df : [Pandas.DataFrame] Data to process 
    ss_vals : [list of float] List of super saturation set points
    param : [str] value to apply the weighted correction to 
    
    Returns
    +++++++
    df: [Pandas.DataFrame] Processed data
    """
    if isinstance(ss_vals,dict):
        df_new =pd.DataFrame()
        for date in list(ss_vals.keys()):
            print(date)
            dt= [pd.to_datetime(d) for d in list(date.split(' - '))]
            data = df[dt[0]: dt[1]]
            print(data)
            print(data.columns.to_numpy())
            x = [ss for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            X_new = [float(list(ss.split("calc_setpt"))[-1])*np.ones(len(data[x[0]].to_numpy())) for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            X_new = np.asarray(X_new, dtype=float)
            X_new = X_new.T
            y = [val for val in data.columns.to_numpy() if f'{param}_avg_setpt' in val]
            y_new = [f'{param}_cor_setpt{list(ss.split("calc_setpt"))[-1]}' for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            X = data[x].to_numpy()
            Y = data[y].to_numpy() 
            X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
            X_new = np.nan_to_num(X_new, nan=0.0, posinf=1e6, neginf=-1e6)
            Y = np.nan_to_num(Y, nan=0.0, posinf=1e6, neginf=-1e6)
            X_scaled = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
            data[y_new] = rowwise_linfit(X,Y,X_new)
            if df_new.empty:
                df_new = data
            else:
                df_new = pd.concat([df_new,data])
        return df_new
    else: 
        x = [ss for ss in df.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
        X_new = [float(list(ss.split("calc_setpt"))[-1])*np.ones(len(df[x[0]].to_numpy())) for ss in df.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
        X_new = np.asarray(X_new, dtype=float)
        X_new = X_new.T
        y = [val for val in df.columns.to_numpy() if f'{param}_avg_setpt' in val]
        y_new = [f'{param}_cor_setpt{list(ss.split("calc_setpt"))[-1]}' for ss in df.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
        X = df[x].to_numpy()
        Y = df[y].to_numpy()
        X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
        X_new = np.nan_to_num(X_new, nan=0.0, posinf=1e6, neginf=-1e6)
        Y = np.nan_to_num(Y, nan=0.0, posinf=1e6, neginf=-1e6)
        X_scaled = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
        df[y_new] = rowwise_linfit(X,Y,X_new)
        return df

"""=====V Correct to STP V====="""
def stp_corr(df,cols,Tstp = 273.15, Pstp= 1013.25):
    r"""
    Applies a STP correction to the number concentration.
    ----------
    Paramaters
    ++++++++++
    df : [Pandas.DataFrame] Data to process.  
    cols : [list of str] Columns to apply correction to.
    Tstd : [float] Standard temperature in K (default = 273.15).
    Pstd : [float] Standard Pressure in hPa(default: 1013.25).

    Returns
    +++++++
    df : [Pandas.DataFrame] Original Dataframe with STP corr columns added
    """
    Tact = (df['T(C)_sample'].to_numpy() + 273.15) #temp in kelvin
    Pact = (df['P(hPA)_sample'].to_numpy())#pressure in hPa
    df['T(C)_stp'] = Tstp-273.15
    df['P(hPA)_stp'] = Pstp
    for col in cols:
        new_col= 'cor_stp'.join(col.split('cor'))
        df[new_col] = df[col]*Pstp/Pact*Tact/Tstp
    return df


def CCN_EBAS(file_in,folder_out, ss_vals):
    """
    Takes in a path to processed CCN file and generates a NASA AMES formated file
    ----------
    Paramaters
    ++++++++++
    file_in : [str/path-like] Path to the processed CCN file
    folder_out : [str/path-like] Path to folder to place EBAS file
    ss_vals: [list of floats] ss% set points

    Returns
    ++++++++++
    NONE
    """
    df = pd.read_csv(file_in)#read in csv skipping first row of verbose column headings
    df=df.set_index('Datetime UTC')

    df = df.fillna(0)
    df.index = pd.to_datetime(df.index)
    dates= df.index.to_list()
    # input(np.isnan(np.sum(df.values.tolist())))
    ccn_corr_cols = [f'N(cm-3)_cor_setpt{ss}' for ss in ss_vals]
    ccn_cols = [f'cloud_condensation_nuclei_number_concentration, 1/cm3, SS={sp}%' for sp in ss_vals]
    data = df[ccn_corr_cols].values.tolist()
    # df['flag'] = [000]
    # df['flag']['Q_flag' ==1] = [662]
    # df['flag']['integrity_flag' == 1] = [111]
    flags = [[[000] for j in range(len(data[i]))] for i in range(len(dates))]
    data =  list(map(list, zip(*data)))
    flags = list(map(list, zip(*flags)))

    date_list = []
    for i in range(len(dates)-1):
        if i == 0:
            date_list.append(pd.to_datetime(dates[i]))
        elif (pd.to_datetime(dates[i+1])-pd.to_datetime(dates[i]) > pd.Timedelta(1, 'hr')):
            date_list.append(pd.to_datetime(dates[i]))
            date_list.append(pd.to_datetime(dates[i+1]))
    date_list.append(pd.to_datetime(dates[-1]))  
    # data = np.nan_to_num(data, nan=1, posinf=1e6, neginf=-1e6)
    # input(np.isnan(np.sum(data)))
    ebas_genfile(folder_out, data, flags, dates, ccn_cols)
