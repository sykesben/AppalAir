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
import re

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

def smps_means(files,freq='d'):
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
        file=file.set_index("DateTime Sample Start") #Set index
        smps = file if smps.empty else pd.concat([smps, file])
    smps.index = pd.to_datetime(smps.index, format='mixed')
    smps = smps[smps['Total Concentration (#/cm³)'].notna()]
    cols = ['Median (nm)',"Mean (nm)",'Geo. Mean (nm)','Mode (nm)','Geo. Std. Dev',]
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
    logdp = np.log10(dp)
    dlogdp = np.diff(logdp)
    dlogdp = np.append(dlogdp, dlogdp[-1])  # pad last bin
    weighted = smps * dlogdp
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

def smps_data_corr(files, freq='D'):
    """
    SMPS processing with:
    - cumulative distribution
    - approximate multiply-charged correction (1+, 2+, 3+)
    - CCN-closure-ready output
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
    """
    smps_total = pd.DataFrame()

    for f in files:
        smps = pd.read_csv(f)
        smps = smps.set_index("DateTime Sample Start")
        smps.index = pd.to_datetime(smps.index, format='mixed')
        smps = smps[smps['Total Concentration (#/cm³)'].notna()]
        numsmps = [s for s in smps.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
        total = smps['Total Concentration (#/cm³)'].to_numpy()
        # numerical sort
        numsmps = sorted(numsmps, key=lambda x: float(x))
        # plt.ion()
        fig, ax = plt.subplots()
        smps = smps[numsmps]
        PSD= smps.copy().mean()
        PSD.index =(PSD.index.to_numpy().astype(float))

        dp = np.array([float(n) for n in numsmps])  # nm
        min_dp = 15  # nm
        dpmsk = dp >= min_dp

        dp = dp[dpmsk]
        numsmps = np.array(numsmps)[dpmsk]

        # update PSD dataframe
        N = smps[numsmps].copy().astype(float)

        logdp = np.log10(dp)
        dlogdp = np.diff(logdp)
        dlogdp = np.append(dlogdp, dlogdp[-1])

        # MULTI-CHARGE CORRECTION (APPROX)
        charge_shift = {2: 1.7,3: 2.4}
        charge_frac = {2: 0.10,3: 0.02}
        corrected = pd.DataFrame(index=smps.index, columns=numsmps,dtype=float)
        corrected.iloc[:, :] = N.values

        for i in range(len(dp) - 1, -1, -1):
            d = dp[i]
            # current corrected concentration
            Ncorr = corrected.iloc[:, i].copy()
            # interpolation of CURRENT corrected PSD
            interp_func = interp1d(dp,corrected.values,axis=1,bounds_error=False,fill_value=0.0)
            subtraction = np.zeros(len(corrected.index))

            for z in [2, 3]:
                dz = d * charge_shift[z]
                # skip if outside scan range
                if dz > dp.max():
                    continue
                # interpolated larger-particle contribution
                Nz = interp_func(dz)
                subtraction += charge_frac[z] * Nz
            corrected.iloc[:, i] = Ncorr - subtraction

        # prevent negatives
        corrected = corrected.clip(lower=0)
        PSD_corr = corrected.copy().mean()
        PSD_corr.index = (PSD_corr.index.to_numpy().astype(float))
        # ax.plot(PSD_corr,label='Corrected PSD')
        # ax.plot(PSD, label ='PSD')
        # ax.legend()
        # ax.set_ylabel('Counts[#]')
        # ax.set_xlabel('Bin [nm]')
        # ax.set_title('PSD charge adjustment')
        # input('Press enter to exit plot...')
        # plt.ioff()

        cols = []
        dp_arr = dp
        for i, d in enumerate(dp_arr):
            mask = dp_arr >= d
            col = f'>{float(d)}nm'
            cols.append(col)
            smps[col] = (corrected.loc[:, np.array(numsmps)[mask]].multiply(dlogdp[mask], axis=1).sum(axis=1))
        smps['Total Concentration (#/cm³)'] = total
        cols.append('Total Concentration (#/cm³)')
        smps_total = smps if smps_total.empty else pd.concat([smps_total, smps])

    smps_total = smps_total.resample(freq).mean()
    smps_total.index.name = 'Date'
    return smps_total, cols

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
    cols = ['T(C)_inlet','T1(C)','T(C)_sample','T(C)_OPC','T(C)_nafion','Q(lpm)_sample','Q(lpm)_sheath','P(hPA)_sample']
    ss_cols = []
    for c in ccn.columns.to_numpy(): 
        if (f'N(cm-3)_cor_setpt' in c) | (f'ss(%)_calc_setpt' in c) | (f'N(cm-3)_avg_setpt' in c):
            cols.append(c)
        if (f'N(cm-3)_cor_setpt' in c):
            ss_cols.append(c)
    ccn = ccn[cols]
    ccn = ccn.resample(freq).mean()
    ccn.index.names = ['Date']
    return ccn,cols,ss_cols

def comb_files(smps_files,ccn_files, freq = 'd', D50 = False,chem = 0, extra_fout = 0):
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
        smps, diam_cols = smps_data_corr(smps_files,freq)
    else:
        smps, diam_cols = smps_data(smps_files,freq)
    # input(smps)
    data = pd.merge(ccn[ccn_cols],smps[diam_cols],left_index = True, right_index = True)
    print(chem)
    if isinstance(chem,str):
        acsm, spec = master_data(chem)
        data = pd.merge(data, acsm ,left_index = True, right_index = True)
    if extra_fout != 0:
        data.to_csv(extra_fout)
    return data,ss_cols,diam_cols
