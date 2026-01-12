"""
Date: 10/31/2025
Author: Ben Sykes
Purpose: Process through raw CCN data and convert it according to ACTRIS formating standard
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import os
import datetime as dt 
import scipy as sp
from scipy.linalg import lstsq
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from CCN_EBAS_convert import ebas_genfile
from pathlib import Path

"""=====I Read in CCN File and CCN.ini file I====="""
def readin(path):


    data = pd.read_csv(path,skiprows=lambda x:x==1)#read in csv skipping first row of verbose column headings
    data=data.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
    data.index = pd.to_datetime(data.index)
    cols_rename = {'particle number concentration (cm-3)': 'N(cm-3)','inlet temperature (°C)': 'T(C)_inlet','temperature of TEC 1 (°C)':'T1(C)','temperature of TEC 2 (°C)':'T2(C)'
                   ,'temperature of TEC 3 (°C)':'T3(C)','sample temperature (°C)': 'T(C)_sample','OPC temperature (°C)':'T(C)_OPC','nafion temperature (°C)':'T(C)_nafion',
                   'sample flow rate (lpm)': 'Q(lpm)_sample','sheath flow (lpm)':'Q(lpm)_sheath','reported supersaturation from onboard instrument calibration (%)':'ss(%)_setpt',
                   'sample pressure (hPa)':'P(hPA)_sample','temperature gradiant setpoint (°C)': 'TG(C)_setpt'}
    data.rename(columns=cols_rename,inplace=True)
    return data,cols_rename

def readini(path):
    date = ''
    TGslope, TGintercept, TGdum = 0,0,0 #temperature gradient values
    date = ''
    ss_list = {'all': [0.1,0.15,0.25,0.4,0.7]} #[0.1,0.15,0.2,0.25,0.4,0.6,0.7,0.8,1.0]
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
        # elif 'SS Table Settings.SS %' in l: 
        #     ss_list.append(float(list(l.split('='))[-1]))
        # elif 'SS Table 1 0.SS Table Settings.SS %' in l: 
        #     ss_list.append(float(list(l.split('='))[-1]))
    file.close()
    ss_list = np.unique(ss_list)
    return TGdum, TGslope, TGintercept, ss_list,date

"""=====II Calculate Super Saturation from Temp Gradient II====="""
def sup_sat(T1,T2,A=16.01,B=1.03,TG=7):
    '''
    Super-Saturation Calculator from Temperature gradient. Assume linear relationship 
    between ss% and TG (for ss%>0.1%) such that TG = A*ss% + B.
    INPUTS: 
    T1: Temperature at TEC 1 [°C] (variable: 'T1_N12')
    T2: Temperature at TEC 2 [°C] (variable: 'T2_N12')
    A: Slope (variable: 'Temp Gradient Slope')
    B: Intercept (variable: 'Temp Gradient Y-Intercept')
    TG: not used

    Outputs:
    tg = temperature gradient len of T1 [C]
    ss: array of super saturation. Len of T1 [%]   
    '''
    tg=(T2-T1)*2 #Use this value according to actris formating
    ss = (tg - B)/A
    return tg,ss



"""=====III Average data over time and super saturation set point III====="""
def time_avg_ss(df, deltat='1h', ss_vals = [], ssflag = True):
    '''
    Groups values by machine set super saturation allowing for group time averaged values
    INPUTS: 
    df: Dataframe to time average [pd.Dataframe]
    deltat: time period to re-average to [str] (default: '1h', must be in '#n' format)
        possible values: m = minute, h = hour, d = day, M = month, y = year
    ss_vals: possible values of ss from machine. If no value is given look at all recorded set points.
        IF a dictionary, the keys should be dates and the values ss lists to apply seperate ss values over
        different date ranges.[list or dict] (default: [])
    ss_flag: Does the measured ss devate more than 20% from the set point. [bool] (default: True)
        Outputs:
    df:  Time averaged Dataframe [pd.Dataframe]
    '''
    if isinstance(ss_vals, dict):
        df_new =pd.DataFrame()
        for date in list(ss_vals.keys()):
            ss_list = ss_vals[date]
            dt= [pd.to_datetime(d) for d in list(date.split(' - '))]
            data = df[dt[0]: dt[1]]
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
            data_new = data.resample(deltat).mean() #Start at the top of the hour,day, month
            times = data_new.index.to_numpy()
            completeness = []
            for i in range(len(times)):
                ts = times[i]
                tf = ts +pd.Timedelta(dt_num, deltat[1])- pd.Timedelta(1,'min')
                slct = data.loc[ts:tf]
                completeness.append(len(slct)/dt) #how many minutes/vs minutes in an hour
                print(slct['ss(%)_setpt'])
                slct= slct[slct['ss(%)_setpt'] == slct['ss(%)_setpt'].shift()]
                input(slct['ss(%)_setpt'])
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
        return df_new
    else: 
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
            # input(f'{len(slct)}/{dt} = {len(sl
            # print(slct['ss(%)_setpt'])
            # slct1= slct[slct['ss(%)_setpt'] == slct['ss(%)_setpt'].shift()]
            # print('updates')
            # print(slct1['ss(%)_setpt'])
            # slct2 = slct.loc[((slct['ss(%)_setpt'] == ss)& (slct['ss_flag']== 0))]
            # print('updates')
            # input(slct2['ss(%)_setpt'])
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
        return df_new


"""=====IV Apply a weighted linear correction to the CCN particle number to ss% IV====="""
def rowwise_linfit(X,X_new, Y):
    """
    Fits Y_i = a_i * X_i + b_i for each row i.
    
    Parameters
    ----------
    X : array-like, shape (n_rows, n_points)
        Independent variable values.
    Y : array-like, shape (n_rows, n_points)
        Dependent variable values.
    
    Returns
    -------
    slopes : ndarray, shape (n_rows,)
        Fitted slopes for each row.
    intercepts : ndarray, shape (n_rows,)
        Fitted intercepts for each row.
    Y_fit : ndarray, shape (n_rows, n_points)
        Fitted values using the row-wise linear models.
    """
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
    return Y_fit

def weighted_corr(df,ss_vals,param):
    """
    Applies a weighted average to generate a linear fit for N vs ss
    INPUTS:
    df: Data to process (DataFrame)
    ss_vals: List of super saturation set points (list)
    param: value to apply the weighted correction to (str)
    
    OUTPUT:
    df: Processed data(DataFrame)   
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
            # X_new = [print(list(ss.split("calc_setpt"))[-1])*np.ones(len(data[x[0]].to_numpy())) for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            X_new = [float(list(ss.split("calc_setpt"))[-1])*np.ones(len(data[x[0]].to_numpy())) for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            X_new = np.asarray(X_new, dtype=float)
            X_new = X_new.T
            y = [val for val in data.columns.to_numpy() if f'{param}_avg_setpt' in val]
            y_new = [f'{param}_cor_setpt{list(ss.split("calc_setpt"))[-1]}' for ss in data.columns.to_numpy() if 'ss(%)_calc_setpt' in ss]
            # input(y_new)
            X = data[x].to_numpy()
            Y = data[y].to_numpy() 
            #Recondition X for linear fit
            X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
            X_new = np.nan_to_num(X_new, nan=0.0, posinf=1e6, neginf=-1e6)
            Y = np.nan_to_num(Y, nan=0.0, posinf=1e6, neginf=-1e6)
            X_scaled = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
            data[y_new] = rowwise_linfit(X,Y,X_new)
            # input(data[y_new])
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
        # input(y_new)
        X = df[x].to_numpy()
        Y = df[y].to_numpy() 
        #Recondition X for linear fit
        X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
        X_new = np.nan_to_num(X_new, nan=0.0, posinf=1e6, neginf=-1e6)
        Y = np.nan_to_num(Y, nan=0.0, posinf=1e6, neginf=-1e6)
        X_scaled = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
        df[y_new] = rowwise_linfit(X,Y,X_new)
        # input(df[y_new])
        return df

"""=====V Correct to STP V====="""
def stp_corr(df, cols):
    temp_act = (df['T(C)_sample'].to_numpy()+273.15) #temp in kelvin
    press_act = (df['P(hPA)_sample'].to_numpy())#pressure in hPa
    temp_stp = 273.15 #standard temp in kelvin = 0C
    press_stp= 1013.25 #standard pressure in hPa
    df['T(C)_stp'] = temp_stp-273.15
    df['P(hPA)_stp'] = press_stp
    for col in cols:
        new_col= 'cor_stp'.join(col.split('cor'))
        df[new_col] = df[col]*press_stp/press_act*temp_act/temp_stp
    return df


file = r"C:\Users\Benjamin Sykes\Documents\CCN_data\CCN_Clean_2025_1min.csv"
ini_file = r"C:\Users\Benjamin Sykes\Documents\CCN_data\CCN 100.ini"
file_out = r"C:\Users\Benjamin Sykes\Documents\CCN_data\CCN_Processed_2025_1hr.csv"
ebas_out = Path(r"C:\Users\Benjamin Sykes\Documents\CCN_data")
df, dt_dct= readin(file)
TGdum, slope, intercept, ss_list, date = readini(ini_file)
ss_vals = [0.1,0.15,0.25,0.4,0.7]

T1= df['T1(C)'].to_numpy()
T2= df['T2(C)'].to_numpy()
df['TG(C)_calc'],df['ss(%)_calc'] = sup_sat(T1,T2,slope, intercept)

#ss flag calculation for the weighted correction
df['ss_dev'] = 100*np.abs(df['ss(%)_calc'].to_numpy() - df['ss(%)_setpt'].to_numpy())/((df['ss(%)_calc'].to_numpy() + df['ss(%)_setpt'].to_numpy())/2)
df['ss_flag'] = (df['ss_dev'].to_numpy()>20.0).astype(int)
ccn_corr_cols = ['N(cm-3)_cor_setpt0.1','N(cm-3)_cor_setpt0.15','N(cm-3)_cor_setpt0.25','N(cm-3)_cor_setpt0.4','N(cm-3)_cor_setpt0.7']
# bins = [c for c in df.columns.to_numpy() if 'Bin' in c]
df = time_avg_ss(df,ss_vals=ss_vals)
df = weighted_corr(df, ss_vals=ss_vals, param='N(cm-3)')
df = stp_corr(df, ccn_corr_cols)
#Drop unnecesary columns
df.dropna(axis =1, how='all', inplace=True)
df.dropna(axis =0,thresh = 5, inplace = True)
df = df.drop(columns=['T2(C)', 'T3(C)','TG(C)_calc','N(cm-3)','ss(%)_setpt', 'TG(C)_setpt','ss(%)_calc', 'ss_dev'])

#Flags
ss_flag = df.pop('ss_flag')# move all flags to end of dataframe
df['ss_flag'] = ss_flag
df['T1_flag1'] = (df['T1(C)'].to_numpy()<30.0).astype(int) # #T1 less than 30 deg and less than Tinlet 
df['integrity_flag'] = df['avg_complete'].to_numpy()>.75
# df['T1_flag2'] = ((np.nanmean(df['T1(C)'].to_numpy())-5.0)<df['T1(C)'].to_numpy()<(np.nanmean(df['T1(C)'].to_numpy())+5.0)).astype(int) #T1 should stay within a 10ish degree band
# df['Q_flag'] =((np.nanmean(df['Q(lpm)_sample'].to_numpy())-5/100*np.nanmean(df['Q(lpm)_sample'].to_numpy()))<df['Q(lpm)_sample'].to_numpy()<(np.nanmean(df['Q(lpm)_sample'].to_numpy())+5/100*np.nanmean(df['Q(lpm)_sample'].to_numpy()))).astype(int) #sample flow should stay within about 5%
# df['N_flag'] = (df['N(cm-3)'].to_numpy()>5000.0).astype(int)
df['date_run'] = pd.to_datetime('now',utc=True).date()
df['date_param'] = pd.to_datetime(date, utc=True).date()
df['ss_slope'] = slope
df['ss_int'] = intercept
df.to_csv(file_out)
dates= df.index.to_list()
ccn_corr_cols = ['N(cm-3)_cor_setpt0.1','N(cm-3)_cor_setpt0.15','N(cm-3)_cor_setpt0.25','N(cm-3)_cor_setpt0.4','N(cm-3)_cor_setpt0.7']
ccn_cols = [f'ccnc[ss={sp}]' for sp in ss_vals]
# input(ccn_cols)
data = df[ccn_corr_cols].values.tolist()
flags = [[[000] for j in range(len(data[i]))] for i in range(len(dates))]
data =  list(map(list, zip(*data)))
flags = list(map(list, zip(*flags)))
# print(data)
# input(np.shape(data))

date_list = []
for i in range(len(dates)-1):
    if i == 0:
        date_list.append(pd.to_datetime(dates[i]))
    elif (pd.to_datetime(dates[i+1])-pd.to_datetime(dates[i]) > pd.Timedelta(1, 'hr')):
        date_list.append(pd.to_datetime(dates[i]))
        date_list.append(pd.to_datetime(dates[i+1]))
date_list.append(pd.to_datetime(dates[-1]))  
# input(date_list)

ebas_genfile(ebas_out, data, flags, dates, ccn_cols)


