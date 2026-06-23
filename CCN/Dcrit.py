"""
Date: 6/23/2026
Author: Ben Sykes
Purpose: Useful plots for calculating critical diameter and kappa values utilizing the CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from CombineData import comb_files,smps_means
from scipy.special import erf
pd.set_option('mode.chained_assignment', None)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.titlesize'] =18

'''======General Functions======+++
General functions used for processing: 
seasons : Produces season name due to corresponding month number
critical_diameter : Calculates a critical diameter from an ss and
  kappa value
kappa_calc : Calculates a hygroscopicity from an ss and Dcrit value
find_mid : Calculates a Critical diameter between two SMPS bin sizes
+++======General Functions======'''
def seasons(num):
    """
    Returns the name of the season a month is part of
    ----------
    Parameters
    ++++++++++
    num : [int] month number

    Returns
    ++++++++++
    season : [str] season name
    """
    if num in [12,1,2]: return 'winter'
    elif num in [3,4,5]: return 'spring'
    elif num in [6,7,8]: return 'summer'
    else: return 'autumn' 

def critical_diameter(ss, kappa=0.1, T=298):
    """
    Estimated critical diameter in nm from Super saturation and given kappa value
    ----------
    Parameters
    ++++++++++
    ss : [float] Supersaturation value [%]
    kappa : [float] hygroscopicity of CCN activation aerosols [1] (default = 0.1)
    T : [float] Temperature of room [K] (default = 298K<-25C)

    Returns
    ++++++++++
    Dcrit : [float] Critical Diameter from CCN/SMPS Calculation [nm]
    """
    sigma = 0.072  # surface tension (N/m)
    Mw = 0.018     # kg/mol
    R = 8.314
    rho_w = 1000   # kg/m3

    A = (4 * sigma * Mw) / (R * T * rho_w)
    ss = float(ss) / 100  # % to fraction
    Dcrit = ((4 * A**3) / (27 * kappa * ss**2))**(1/3)
    return Dcrit * 1e9  # m to nm

def kappa_calc(ss, Dcrit, T=298):
    """
    Calculated kappa from Supersaturation and Critical Diameter
    ----------
    Parameters
    ++++++++++
    ss : [float] Supersaturation value [%]
    Dcrit : [float] Critical Diameter from CCN/SMPS Calculation [nm]
    T : [float] Temperature of room [K] (default value =298K<-25C)

    Returns
    ++++++++++
    kappa : [float] hygroscopicity of CCN activation aerosols [1]
    """
    sigma = 0.072  # surface tension (J/m^2)
    Mw = 0.018     # kg/mol
    R = 8.314      # J/Kmol
    rho_w = 1000   # kg/m3

    A = (4 * sigma * Mw) / (R * T * rho_w)
    ss = float(ss) / 100  # % to fraction
    Dcrit = Dcrit / 1e9   # nm to m
    kappa = (4 * A**3) / (27 * Dcrit**3 *np.log(1+ss)**2)
    kappa = np.where(Dcrit!=0, kappa, 100*np.ones(len(kappa)))
    return kappa

def find_mid(D_top, D_bot, diff_top, diff_bot):
    """
    Calculates the "actual" diameter by finding the weighted midpoint between
    the two nearest diameter bins from the SMPS. 
    ----------
    Parameters
    ++++++++++
    D_top : [float] upper diameter bin size for Critical Diameter [nm]
    D_bot : [float] lower diameter bin size for Critical Diameter [nm]
    diff_top : [float] difference between N_CCN and N_Dtop
    diff_bot : [float] difference between N_CCN and N_Dbot
    
    Returns
    ++++++++++
    diameter : [float] critical diameter [nm]
    """
    delta_D = abs(D_top-D_bot)
    delta_diff = abs(diff_bot) + abs(diff_top)
    dev_from_bottom = abs(diff_bot)/delta_diff
    diameter = delta_D*dev_from_bottom + D_bot
    return diameter

'''======Finding Calculation Functions======+++
The algorhythmic functions used for processing datasets: 
find_cutoff : calculates the critical diameter and hygroscopicity
+++======Finding Calculation Functions======'''
'''Sigmoid function fitting'''
def activation(D50, D, sigma=32):
    """
    Fits an activation curve for the PSD
    ----------
    Parameters
    ++++++++++
    D50 : [float] Critical diameter
    D : [float] Particle size distribution 
    sigma : [float] spread of activation thresholds (default = 32)

    Returns
    ++++++++++
    A : [float] activation curve 
    """
    return 0.5 * (1 + erf((D-D50)/(np.sqrt(2)*sigma)))

'''Find D50'''
def find_Diameter(data, diams, ss_cols,sigma = 32):
    """
    Calculates the critical diameter and hygroscopicity using N_CCN and N_CN data from the SMPS
    ----------
    Parameters
    ++++++++++
    data : [panda Dataframe] combined datasets from CCN and SMPS
    diams : [list of str] diameter bin sizes from SMPS 
    ss_cols : [list of str] N_(@ss) column names from CCN

    Returns
    ++++++++++
    cut_off : [dict] Dictionary of Critical Diameters with keys = ss% // values = [D_crit, Dev_bottom, Dev_top]
    kappa : [dict] Dictionary of hygroscopicity with keys = ss% // values = [kappa, Dev_bottom, Dev_top]
    """
    Dcrits = {}
    kappas = {}

    # Convert to diameters 
    dp = np.array([float(n) for n in diams])
    min_dp = 15  # nm
    dpmsk = dp >= min_dp
    dp = dp[dpmsk]
    logdp = np.log10(dp)
    dlogdp = np.diff(logdp)
    dlogdp = np.append(dlogdp, dlogdp[-1])
    diams = np.array(diams)[dpmsk]          #shortened SMPS columns
    N = data[diams].copy().astype(float)    #pull out just the SMPS columns from data
    # Generate a grf possible D50s 
    D50_grid = np.arange(20, 300, 0.5)
    # Using the possible D50s calculate a grid from the fitted sigmoid function
    A  = activation(D50_grid[None,:],dp[:,None],sigma) # activation grid
    dist = N.values # SMSP values (dt) x (dD) 
    # Using the distribution and activation grid and Dlog(dp) values calculate and estimated CCN distribution grid 
    N_CCN_pred = (dist[:,:,None]*A[None,:,:]*dlogdp[None,:,None]).sum(axis=1)
    # itterate through each ss set point to pull out N_CCN(ss)
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        N_CCN_obs = data[col].to_numpy()  # N_CCN(ss) (dt) x (dD) 
        # find the D50 values that minimizes the deviation between the Predicted and Observed CCN values
        D50_idx = np.argmin(np.abs(N_CCN_pred- N_CCN_obs[:,None]),axis=1) 
        # Pull out the D50 from each index
        D50 = D50_grid[D50_idx]
        # Calculate Kappa from each D50
        kappa = kappa_calc(ss,D50)
        Dcrits[ss] = D50
        kappas[ss] = kappa
    return Dcrits, kappas

'''Find Fact'''
def find_Activation(data, ss_cols):
    Fact_dict ={}
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        Fact_dict[f"Fact(ss={ss})"] = data[col].to_numpy()/data['Total Concentration (#/cm³)'].to_numpy()
    return Fact_dict

"""Calls in the functions if this script is run"""
if __name__ == '__main__':
    '''1=======User Inputs=======1+++
    Supplied inputs for the processing
     - dates to remove
     - input files
     - frequency
     - sigma of activation function (suggested 16-32)
    +++1=======User Inputs=======1'''
    # Known dates to remove from processing
    bad_dates = [[pd.to_datetime('06/01/2024 00:00:00'),pd.to_datetime('07/30/2024 00:00:00')],[pd.to_datetime('10/20/2024 00:00:00'),pd.to_datetime('10/27/2024 00:00:00')],
                 [pd.to_datetime('10/01/2025 00:00:00'),pd.to_datetime('12/31/2025 00:00:00')],[pd.to_datetime('07/01/2025 00:00:00'),pd.to_datetime('07/24/2025 00:00:00')],[pd.to_datetime('08/06/2025 00:00:00'),pd.to_datetime('08/14/2025 00:00:00')],
                 [pd.to_datetime('01/01/2026 00:00:00'),pd.to_datetime('01/09/2026 18:00:00')]]
    # Files to combine
    smps =[r"C:\Users\bensy\Documents\Research\2026_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv"]  #list(input('Provide paths to SMPS file(s). Seerate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN\Processed\CCN_lvl2_2026_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN\Processed\CCN_lvl2_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN\Processed\CCN_lvl2_2024_1hr.csv"]   #r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv"  #[r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv"]#  #list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    f = 'h' # set frequency
    sgm = 16 # 16 seems to minimize noise of data, but between 16-32 seems to work well

    '''2=======Read In=======2+++
    Read in files
     - Combine SMPS and CCN files using comb_files function
     - Pull the non-PSD data from the SMPS file and merge 
     - Drop bad data 
    +++2=======Read In=======2'''
    # Combined files and clean
    data,ss_cols,diams = comb_files(smps,ccn,freq=f)
    smps_mean = smps_means(smps,freq = f) # read in SMPS secondary data
    data = pd.merge(data,smps_mean,left_index = True, right_index = True, how='left')
    mask = pd.Series(False, index=data.index) # drop bad dates
    for date in bad_dates:
        mask |= (data.index >= date[0]) & (data.index <= date[-1])
    data = data[~mask]
    data = data.dropna(axis='index',how='all', subset=diams) # drop empty rows
    data = data.dropna(axis='columns',how='all')

    '''3=======Calculations=======3+++
     Use the find_Diameter function to find critical diameter
    and hygroscopicity
     Us the find_Activation function to calculate activation 
    fractions 
    +++3=======Calculations=======3'''
    # Find critical diameters and hygroscopicity
    Dcrits,Kappas = find_Diameter(data,diams, ss_cols,sigma=sgm)
    dates = data.index.to_numpy()
    # Find activation
    Fact_d = find_Activation(data,ss_cols)

    '''4=======Gen Outputs=======4+++
    Generate the outputs from the calculations
     - From Dcrit dictionary, generate a Dc dataframe
        + set a datetime index
        + rename columns such as "Dcrit(ss=#)"
        + remove unphysical Dcrit values
        + calculate the standard deviation post cleaning
     - From the Kappas dictionary,generate a K dataframe
        + set datetime index
        + rename columns such as "Kappa(ss=#)"
        + remove unphysical Kappa values
     - From the Fact_d, generate a Fact datafram
        + set a datetime index
     - Combine the Dc and K DataFrame, use the overlapping
       dates from both DataFrames to generate a larger output.
     - Merge extra data into this larger DataFrame using the 
       dates of the larger DataFrame.
     - Output DataFrame to specified path as a csv
    +++4=======Gen Outputs=======4'''
    # Critical diameter dataframe
    Dc = pd.DataFrame(Dcrits)
    Dc['Datetime'] = dates
    Dc = Dc.set_index("Datetime")
    Dc.index = pd.to_datetime(Dc.index, format='mixed')
    Dc = Dc.add_prefix('Dcrit(ss=', axis='columns')
    Dc = Dc.add_suffix(')', axis='columns')
    # remove data with physically impossible values 
    Dc[Dc < 10] = np.nan
    Dc[Dc > 300] = np.nan
    DC_STD = Dc.std()
    # Hygroscopicity dataframe
    K = pd.DataFrame(Kappas)
    K['Datetime'] = dates
    K = K.set_index("Datetime")
    K.index = pd.to_datetime(K.index, format='mixed')
    K = K.add_prefix('Kappa(ss=', axis='columns')
    K = K.add_suffix(')', axis='columns')
    K[K>1.5] = np.nan
    Dc[K>1.5] = np.nan
    # remove the worst spikes in the data 
    K['Kappa(ss=0.7)'][K['Kappa(ss=0.7)']>0.4] = np.nan
    K['Kappa(ss=0.4)'][K['Kappa(ss=0.4)']>0.6] = np.nan
    K['Kappa(ss=0.25)'][K['Kappa(ss=0.25)']>0.6] = np.nan
    K['Kappa(ss=0.15)'][K['Kappa(ss=0.15)']>0.6] = np.nan
    # generate the activation Fraction Dataframe
    Fact = pd.DataFrame(Fact_d)
    Fact['Datetime'] = dates
    Fact = Fact.set_index('Datetime')
    Fact.index = pd.to_datetime(Fact.index, format='mixed')
    # Generate Data Out
    Data_out = pd.merge(Dc,K,left_index = True, right_index = True) 
    ccn_data = data[['N(cm-3)_cor_setpt0.1','N(cm-3)_cor_setpt0.15','N(cm-3)_cor_setpt0.25',"N(cm-3)_cor_setpt0.4",'N(cm-3)_cor_setpt0.7']].copy()
    Data_out = pd.merge(Data_out,smps_mean,left_index = True, right_index = True,how='left')
    Data_out = pd.merge(Data_out,ccn_data,left_index = True, right_index = True,how='left')
    Data_out = pd.merge(Data_out,Fact,left_index = True, right_index = True,how='left')
    out = input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '' : Data_out.to_csv(out)
    
    print('finished')