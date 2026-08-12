"""
Date: 4/25/2026
Author: Ben Sykes
Purpose: Call comb_files to combine data from different datasets
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr, kstest
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import least_squares as LSfit
from scipy.optimize import curve_fit
from datetime import datetime 
from scipy.signal import find_peaks
from scipy.special import erf
import re

def logNfit(y, Dp,gstd =1.75,med=80):
    y = np.asarray(y,dtype=float)
    mask = (~np.isnan(y)) & (y>0)
    if mask.sum() < 5:
        return pd.Series({"N_Peaks": np.nan,"R2": np.nan, 'N': np.nan,
                          "Mode1_N": np.nan,"Mode1_CMD": np.nan,"Mode1_GSD": np.nan,
                          "RightMode_N": np.nan,"RightMode_CMD": np.nan,"RightMode_GSD": np.nan})
    logDp = np.log10(Dp)
    x = logDp[mask]
    y = y[mask]
    y_smooth = pd.Series(y).rolling(7, center=True).mean().to_numpy() #smooth out y
    y_smooth[np.isnan(y_smooth)] = y[np.isnan(y_smooth)]
    y = y_smooth
    peaks, props = find_peaks(y,prominence=0.15*y.max())
    CMD_guess = Dp[mask][peaks]
    N_guess = y[peaks]
    p0 = [y.max(),med,gstd]
    pars0, cov0 = curve_fit(lognormal,x,y,p0=p0,bounds=([0,5,1.1],[np.inf,1000,3.5]))
    N1, cmd1, gsd1 = pars0
    fit0 = lognormal(x,pars0[0], pars0[1], pars0[2])
    N_calc = np.trapezoid(fit0,x)
    ss_res = np.sum((y-fit0)**2)
    ss_tot = np.sum((y-y.mean())**2)
    r2 = 1 - ss_res/ss_tot
    if len(peaks) >1:
        try:
            p1 = [y[peaks][0] * np.log10(1.4) * np.sqrt(2*np.pi),Dp[mask][peaks][0],1.4, #left mode
                    y[peaks][1] * np.log10(1.3) * np.sqrt(2*np.pi),Dp[mask][peaks][1],1.3] #right mode
            pars2, cov2 = curve_fit(logbinormal, x, y, p0=p1)
            N1, cmd1, gsd1, N2, cmd2, gsd2 = pars2
            bimodal = logbinormal(x, N1, cmd1, gsd1, N2, cmd2, gsd2)
            right= lognormal(x, N2, cmd2, gsd2)
            N_calc = np.trapezoid(right,x)
            ss_res = np.sum((y-bimodal)**2)
            ss_tot = np.sum((y-y.mean())**2)
            r2 = 1 - ss_res/ss_tot
            return pd.Series({"N_Peaks": len(peaks),"R2": r2, 'N':N_calc,
                            "Mode1_N": N1,"Mode1_CMD": cmd1,"Mode1_GSD": gsd1,
                            "RightMode_N": N2,"RightMode_CMD": cmd2,"RightMode_GSD": gsd2})
        except:
            return pd.Series({"N_Peaks": len(peaks),"R2": r2,'N':N_calc,
                                      "Mode1_N": N1,"Mode1_CMD": cmd1,"Mode1_GSD": gsd1,
                                      "RightMode_N": N1,"RightMode_CMD": cmd1,"RightMode_GSD": gsd1})
    else: 
        return pd.Series({"N_Peaks": len(peaks),"R2": r2,'N':N_calc,
                          "Mode1_N": N1,"Mode1_CMD": cmd1,"Mode1_GSD": gsd1,
                          "RightMode_N": N1,"RightMode_CMD": cmd1,"RightMode_GSD": gsd1})

def lognormal(logDp, N, CMD, GSD):
    sigma = np.log10(GSD)
    return (N/(sigma*np.sqrt(2*np.pi))* np.exp(-(logDp-np.log10(CMD))**2/(2*sigma**2)))

def logbinormal(logDp, N1, cmd1, gsd1, N2, cmd2, gsd2):
    mode1 = lognormal(logDp, N1, cmd1, gsd1)
    mode2 = lognormal(logDp, N2, cmd2, gsd2)
    return mode1 + mode2

def AIC(rss,n,k):
    return n*np.log10(rss/n)+2*k

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

def smps_means(files,freq='d', fit = False):
    '''
    Takes in a list of smps files and returns some values
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
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
    if fit:
        smps[['N_Peaks','LogNormFit_R2', 'N',
                "Mode1_N","Mode1_CMD","Mode1_GSD",
                "RightMode_N","RightMode_CMD","RightMode_GSD"]]= smps.copy().apply(lambda x: logNfit(x[numsmps].values,Dp,x['Geo. Std. Dev'],x['Median (nm)']),axis='columns', result_type="expand")
        smps['Total Concentration Fit(#/cm³)'] = smps['N']
    smps = smps[smps['Total Concentration (#/cm³)'].notna()]
    cols = ['Median (nm)',"Mean (nm)",'Geo. Mean (nm)','Mode (nm)','Geo. Std. Dev','Total Concentration (#/cm³)']
    if fit:
        cols.extend(['Total Concentration Fit(#/cm³)','N_Peaks','LogNormFit_R2','N',"Mode1_N","Mode1_CMD","Mode1_GSD","RightMode_N","RightMode_CMD","RightMode_GSD"])
    smps = smps[cols]
    return smps

def smps_data(files,freq='d'):
    '''
    Takes in a list of smps files and filters out the particle size concentration depending
    on comparable ss% values.
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    smps = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        file=file.set_index("DateTime Sample Start") #Set index
        smps = file if smps.empty else pd.concat([smps, file])
    smps.index = pd.to_datetime(smps.index, format='mixed')

    numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]
    total = smps['Total Concentration (#/cm³)'].to_numpy()
    # numerical sort
    numsmps = sorted(numsmps, key=lambda x: float(x))
    smps = smps[numsmps]

    # Convert to diameters 
    dp = np.array([float(n) for n in numsmps])
    min_dp = 15  # nm
    dpmsk = dp >= min_dp

    dp = dp[dpmsk]
    numsmps = np.array(numsmps)[dpmsk]
    N = smps[numsmps].copy().astype(float)

    logdp = np.log10(dp)
    dlogdp = np.diff(logdp)
    dlogdp = np.append(dlogdp, dlogdp[-1])  # pad last bin
    weighted = N * dlogdp
    cols = []
    for n in numsmps:
        col = f'>{float(n)}nm'
        cols.append(col)
        # bins greater than threshold
        mask = dp >= float(n)
        # apply weighted sum instead of raw sum
        smps[col] = weighted.loc[:, np.array(numsmps)[mask]].sum(axis=1)
    smps['Total Concentration (#/cm³)'] = total
    cols.append('Total Concentration (#/cm³)')
    smps = smps.resample(freq).mean()
    smps.index.names = ['Date']
    return smps,cols

def AE33_data(files,freq = 'h'):
    ae33 = pd.DataFrame()
    for i in range(len(files)): #read in ae33 files and combine
        f = files[i]
        file =pd.read_csv(f) #read in ae33 file
        file=file.set_index("DateTimeUTC") #Set index
        ae33 = file if ae33.empty else pd.concat([ae33, file])
    ae33.index = pd.to_datetime(ae33.index, format='mixed')
    ae33.index.rename('Datetime(UTC)', inplace=True)
    ae33 = ae33.resample(freq).mean()
    # input(ae33)
    return ae33

def smps_data_corr(files, freq='D', fit =False):
    """
    Clean and output the PSD produced by the SMPS

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame
    fit : [str] whether to fit PSD to lognormal dist, and whether to 
                fit only right mode, or full fit (defalt = False(no fit))
                [options: False, 'right', 'full']

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    """
    smps_out = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        smps = pd.read_csv(f) #read in smps file\
        print(smps.columns)
        try:
            smps = smps.set_index("DateTime Sample Start") #Set index
        except:
            smps = smps.set_index("Datetime(UTC)") #Set index
        smps.index = pd.to_datetime(smps.index, format='mixed')
        numsmps = [s for s in smps.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
        numsmps = sorted(numsmps, key=lambda x: float(x))
        nums = smps[numsmps].copy()
        Dp = nums.columns.astype(float).to_numpy()
        if (fit == 'right')|(fit == 'full'):
            smps[['N_Peaks','LogNormFit_R2','N',
                    "Mode1_N","Mode1_CMD","Mode1_GSD",
                    "RightMode_N","RightMode_CMD","RightMode_GSD"]]= smps.copy().apply(lambda x: logNfit(x[numsmps].values,Dp,x['Geo. Std. Dev'],x['Median (nm)']),axis='columns', result_type="expand")
            fits = smps[['N_Peaks','LogNormFit_R2',"Mode1_N","Mode1_CMD","Mode1_GSD","RightMode_N","RightMode_CMD","RightMode_GSD"]]
        # numerical sort
        smps = smps[numsmps]
        new_nums = {}
        for num in numsmps:
            try:
                # ensure all bin diameters are 2 decimal places to correctly merge
                float_num = float(num)
                new_nums[num] = f"{float_num:.2f}"
            except:
                new_nums[num] = num
        smps = smps.rename(columns=new_nums)
        if fit =='right':
            fit_smps = pd.DataFrame(index=smps.index,columns=list(new_nums.values()),dtype=float)
            for date, row in fits.iterrows():
                if np.isnan(row["RightMode_N"]):
                    print((row["RightMode_N"]))
                    continue
                fit_smps.loc[date] = np.log10(10) * lognormal(np.log10(Dp),row["RightMode_N"],row["RightMode_CMD"],row["RightMode_GSD"])
            smps = fit_smps
        #Remove the duplicate headers
        smps = smps.T.groupby(level=0).mean().T
        smps = smps.resample(freq).mean() 
        smps.index.names = ['Date']
        smps_out = smps if smps_out.empty else pd.concat([smps_out, smps])
    numsmps = np.asarray(list(new_nums.values()))
    return smps_out, numsmps

def ccn_data(files, freq ='d'):
    '''
    Takes in a list of CCN files and returns a processed dataframe with important columns 
    for plotting or further analysis
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrame (default = 'd')

    Returns
    ++++++++++
    ccn : [DataFrame] Combined CCN data from all inputted files
    cols : [list of str] Names of used columns from CCN output
    '''
    ccn = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in ccn file
        try:
            file=file.set_index('Datetime(UTC)') #Set index
        except:
            try:
                file=file.set_index('Datetime UTC') #Set index
            except:
                file = file.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
        file.index = file.index.rename('Datetime(UTC)')
        if i == 0:
            ccn = file
        else:
            ccn = pd.concat([ccn,file])
    ccn.index = pd.to_datetime(ccn.index, format='mixed')
    cols = ['T(C)_inlet','T1(C)','T(C)_sample','T(C)_OPC','T(C)_nafion','Q(lpm)_sample','Q(lpm)_sheath','P(hPa)_sample']
    ss_cols = []
    for c in ccn.columns.to_numpy(): 
        if (f'N(cm-3)_cor_setpt' in c) | (f'ss(%)_calc_setpt' in c) | (f'N(cm-3)_avg_setpt' in c):
            cols.append(c)
        if (f'N(cm-3)_cor_setpt' in c):
            ss_cols.append(c)
    print(ccn.columns)
    ccn = ccn[cols]
    ccn = ccn.resample(freq).mean()
    ccn.index.names = ['Date']
    return ccn,cols,ss_cols

def comb_files(smps_files,ccn_files, freq = 'd', D50 = True,chem = 0, extra_fout = 0,smps_fit = False):
    '''
    Takes in a list of CCN and SMPS files and returns a combined dataframe with important columns 
    from both for plotting or further analysis
    ----------

    Parameters
    ++++++++++
    smps_files : [list of str] Paths to SMPS files
    ccn_files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrames (default = 'd')
    ss : [list of floats] ss% set points from CCN (default = [0.1,0.7])
    chem : [list of str] Paths to Chemistry data if included (default = 0)
    cortype : [float] if specific correction column is wanted (default = '')

    Returns
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data from all inputted files
    '''
    ccn, ccn_cols,ss_cols = ccn_data(ccn_files,freq)
    if D50: 
        smps, diam_cols = smps_data_corr(smps_files,freq,fit=smps_fit)
    else:
        smps, diam_cols = smps_data(smps_files,freq)
    data = pd.merge(ccn[ccn_cols],smps[diam_cols],left_index = True, right_index = True)
    if isinstance(chem,str):
        acsm, spec = master_data(chem)
        data = pd.merge(data, acsm ,left_index = True, right_index = True)
    if extra_fout != 0:
        data.to_csv(extra_fout)
    return data,ss_cols,diam_cols

# def MAF_calc(vol_files,ae33_files,freq = 'd'):
#     bc = AE33_data(ae33_files,freq)
#     vol = 