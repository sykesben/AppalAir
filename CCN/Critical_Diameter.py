"""
Date: 4/25/2026
Author: Ben Sykes
Purpose: Useful plots for calculating critical diameter and kappa values utilizing the CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import least_squares as LSfit
from datetime import datetime 
from CombineData import comb_files
import re
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
    sigma = 0.072  # surface tension (N/m)
    Mw = 0.018     # kg/mol
    R = 8.314
    rho_w = 1000   # kg/m3

    A = (4 * sigma * Mw) / (R * T * rho_w)
    ss = float(ss) / 100  # % to fraction
    Dcrit = Dcrit / 1e9   # nm to m
    if Dcrit ==0.0:
        kappa = 100
    else:
        kappa = (4 * A**3) / (27 * Dcrit**3 *ss**2)
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
find_perc_above : Calculates the fraction of particles that are
  larger than Dcrit and have activated as CCN
find_activation : Calculates the activation fraction of aerosols
  larger than a given diameter
+++======Finding Calculation Functions======'''
def find_cutoff(data, diams, ss_cols, Kappa = 0.1):
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
    cut_off = {}
    kappa = {}
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        cnt_at_ss = data[col]  # N_CCN(ss) 
        devs = []
        for diam in diams:
            cnt_at_diam = data[diam] # N_CN
            dev = cnt_at_ss-cnt_at_diam # find where 100% of the total particles can activate as CCN
            devs.append(float(dev))
        devs = np.array(devs)
        devs = devs[:-1] #remove the 'total particle' value from end of array
        indices = np.where(~np.isnan(devs))[0]  # Get integer indices
        devs = devs[indices]  # Use integer indices to filter 

        notFound = True
        start = 0
        end = len(devs)-1
        while notFound: # Binary sort through deviations to find SMPS bin size
            mid = (start + end)//2 
            s_dev = devs[start]
            e_dev = devs[end]
            m_dev = devs[mid]
            if s_dev >0: #If the critical diameter is the lowest diameter measured from the SMPS
                notFound = False
                diam_top = float(diams[start].replace('>','').replace('nm','')) 
                kappa_top = kappa_calc(ss,diam_top)
                diam_bot = 0
                kappa_bot = kappa_calc(ss,diam_bot)
                diam_mid = find_mid(diam_top,diam_bot, float(devs[start]),float(devs[start]))
                kappa_mid = kappa_calc(ss,diam_mid)
                cut_off[ss] = [diam_mid, diam_bot,diam_top]
                kappa[ss] = [kappa_mid, kappa_bot,kappa_top]
            if  m_dev * s_dev > 0: # if the middle and start are the same sign
                start = mid
            elif m_dev * e_dev > 0: # if the middle and end are the same sign
                end = mid
            if abs(start-end) == 1: # Once the top and bottom bin sizes are adjacent
                notFound = False
                diam_top = float(diams[end].replace('>','').replace('nm',''))
                kappa_top = kappa_calc(ss,diam_top)
                diam_bot = float(diams[start].replace('>','').replace('nm',''))
                kappa_bot = kappa_calc(ss,diam_bot)
                diam_mid = find_mid(diam_top,diam_bot, float(devs[end]),float(devs[start]))
                kappa_mid = kappa_calc(ss,diam_mid)
                cut_off[ss] = [diam_mid, diam_bot,diam_top]
                kappa[ss] = [kappa_mid, kappa_bot,kappa_top]
    return cut_off, kappa

def find_perc_above(data, diams, ss_cols, grtr_dict, Kappa = 0.1):
    """
    Calculates the number of particles larger than a kohler calculated critical diameter in comparison to the
    CCN activation fraction
    ----------
    Parameters
    ++++++++++
    data : [panda Dataframe] combined datasets from CCN and SMPS
    diams : [list of str] diameter bin sizes from SMPS 
    ss_cols : [list of str] N_(@ss) column names from CCN
    grtr_dict : [dict] Dictionary containing Date/Fraction of Aerosols larger than the Critical Diameter/Activation Fraction
    Kappa : [float] hygroscopicity value to compare against [1]

    Returns
    ++++++++++
    grtr_dict : [dict] Dictionary containing Dates/Fraction of Aerosols larger than the Critical Diameter/Activation Fraction
    """
    dt = pd.to_datetime(data.name)
    year = str(dt.year).replace('20','')
    month = dt.month
    grtr_dict['Date'].append(f'{month}/{year}')
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        est_D = critical_diameter(ss,kappa=Kappa) # Kohler defined critical Diameter at given hygoscopicity and super saturation
        greater_than_diam = [d for d in diams[:-1] if float(d.replace('>','').replace('nm',''))>est_D][0] # first bin size greater than calculated Diameter
        frac_above_diam = data[greater_than_diam]/data[diams[-1]] # N_>D/N_CN
        frac_at_ss = data[col]/data[diams[-1]] # N_CCN/N_CN
        grtr_dict[f'F>Dcrit({ss})'].append(frac_above_diam)
        grtr_dict[f'Fact({ss})'].append(frac_at_ss)
    return grtr_dict

def find_activation(data, diams, ss_cols):
    """
    Calculates the activation fraction of aerosols larger than a given diameter
    ----------
    Parameters
    ++++++++++
    data : [panda Dataframe] combined datasets from CCN and SMPS
    diams : [list of str] diameter bin sizes from SMPS 
    ss_cols : [list of str] N_(@ss) column names from CCN

    Returns
    ++++++++++
    act_perc : [dict] Dictionary containing the activation fraction at each diameter
    """
    act_perc = {}
    d_index = [float(x.replace('>','').replace('nm','')) for x in diams[:-1]]
    act_perc['Diameter'] = d_index
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        cnt_at_ss = data[col]
        Acts = []
        for diam in diams:
            cnt_at_diam = data[diam]
            Act = cnt_at_ss/cnt_at_diam #% of the total particles that activated as CCN
            Acts.append(float(Act))
        Acts = np.array(Acts)
        Acts = Acts[:-1]
        indices = np.where(~np.isnan(Acts))[0]  # Get integer indices
        Acts = Acts[indices]  # Use integer indices to filter 
        act_perc[f'F_act_{ss}'] = Acts 
    return act_perc

'''======Generate Curves From the Data Sets======+++
Functions used for generating plots at different frequencies:
cut_off_curve : Generates critical diameter curves
kappa_curve : Generates hygroscopicity curves
Fact_curve : Activation fraction at diameter curves
Fract_curve : Activation fraction and fraction above Dcrit
+++======Generate Curves From the Data Sets======'''
def cut_off_curve(data, ideal = 0,freq = 'M'):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots of the critical diameter
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    ideal : [list of str] If idealized columns included in dataset, pass a list of ideal column names (default = 0 = no value)
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
        -> if the time resolution of the input data is lower than monthly, set freq to 'M'

    Returns
    ++++++++++
    none 
    '''
    cols = [c for c in data.columns.to_numpy() if 'ideal' not in c] 
    if freq =='S':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        month_vals = np.array([(c.split('_')[-1]).split('/')[0] for c in cols])
        szns= [seasons(float(c)) for c in month_vals]
        Scols = [' - '.join([seasons(float(x)), y]) for x, y in zip(month_vals,year_vals)]
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        if ideal != 0:
            szn_cols = ideal
        else:
            szn_cols = []
        for yr in np.unique(year_vals):
            for szn in np.unique(szns):
                s_idx = [i for i,c in enumerate(Scols) if (szn in c)&(yr in c)]
                cut_col = [cols[i] for i in list(set(s_idx)&set(cutoffs))]
                low_col = [cols[i] for i in list(set(s_idx)&set(lowers))]
                up_col = [cols[i] for i in list(set(s_idx)&set(uppers))]  
                data[f'D_cutoff_{szn}/{yr}'] = data[cut_col].mean(axis =1) # calculate a seasonal mean
                data[f'D_lower_{szn}/{yr}'] = data[low_col].mean(axis =1)
                data[f'D_upper_{szn}/{yr}'] = data[up_col].mean(axis =1)
                szn_cols.extend([f'D_cutoff_{szn}/{yr}',f'D_lower_{szn}/{yr}',f'D_upper_{szn}/{yr}'])
        data = data[szn_cols]
    if freq =='Y':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        if ideal != 0:
            yr_cols = ideal
        else:
            yr_cols = []
        for yr in np.unique(year_vals):
            y_idx = [i for i,c in enumerate(year_vals) if yr in c]
            cut_col = [cols[i] for i in list(set(y_idx)&set(cutoffs))]
            low_col = [cols[i] for i in list(set(y_idx)&set(lowers))]
            up_col = [cols[i] for i in list(set(y_idx)&set(uppers))]            
            data[f'D_cutoff_{yr}'] = data[cut_col].mean(axis =1) # caclulate a yearly mean
            data[f'D_lower_{yr}'] = data[low_col].mean(axis =1)
            data[f'D_upper_{yr}'] = data[up_col].mean(axis =1)
            yr_cols.extend([f'D_cutoff_{yr}',f'D_lower_{yr}',f'D_upper_{yr}'])
        data = data[yr_cols]
    # months and seasons to remove
    bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
    bad_mnths = ['9/25', '10/25','11/25']
    slct = [c for c in data.columns.to_numpy() if ('cutoff' in c)&(not any(y in c for y in bad_szns))&(not any(y in c for y in bad_mnths))]# & ('autumn' not in c)] 
    scat_call(data, slct, freq=freq)
    input(data)

def kappa_curve(data,freq = 'M', plot = ['scat']):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots of the hygroscopicity
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
        -> if the time resolution of the input data is lower than monthly, set freq to 'M'
    plot : [list of str] type of plot to generate (default = ['scat'])
        + 'scat' - scatter plot
        + 'line' - line plot

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    if freq =='S':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        month_vals = np.array([(c.split('_')[-1]).split('/')[0] for c in cols])
        szns= [seasons(float(c)) for c in month_vals]
        Scols = [' - '.join([seasons(float(x)), y]) for x, y in zip(month_vals,year_vals)]
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        szn_cols = []
        for yr in np.unique(year_vals):
            for szn in np.unique(szns):
                s_idx = [i for i,c in enumerate(Scols) if (szn in c)&(yr in c)]
                cut_col = [cols[i] for i in list(set(s_idx)&set(cutoffs))]
                low_col = [cols[i] for i in list(set(s_idx)&set(lowers))]
                up_col = [cols[i] for i in list(set(s_idx)&set(uppers))]    
                data[f'Kappa_cutoff_{szn}/{yr}'] = data[cut_col].mean(axis =1)
                data[f'Kappa_lower_{szn}/{yr}'] = data[low_col].mean(axis =1)
                data[f'Kappa_upper_{szn}/{yr}'] = data[up_col].mean(axis =1)
                szn_cols.extend([f'Kappa_cutoff_{szn}/{yr}',f'Kappa_lower_{szn}/{yr}',f'Kappa_upper_{szn}/{yr}'])
        data = data[szn_cols]
    if freq =='Y':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        yr_cols = []
        for yr in np.unique(year_vals):
            y_idx = [i for i,c in enumerate(year_vals) if yr in c]
            cut_col = [cols[i] for i in list(set(y_idx)&set(cutoffs))]
            low_col = [cols[i] for i in list(set(y_idx)&set(lowers))]
            up_col = [cols[i] for i in list(set(y_idx)&set(uppers))]            
            data[f'Kappa_cutoff_{yr}'] = data[cut_col].mean(axis =1)
            data[f'Kappa_lower_{yr}'] = data[low_col].mean(axis =1)
            data[f'Kappa__upper_{yr}'] = data[up_col].mean(axis =1)
            yr_cols.extend([f'Kappa_cutoff_{yr}',f'Kappa_lower_{yr}',f'Kappa_upper_{yr}'])
        data = data[yr_cols]
    bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
    slct = [c for c in data.columns.to_numpy() if ('cutoff' in c)&(not any(y in c for y in bad_szns))]
    if (freq !='d')&('scat' in plot):
        scat_call(data, slct, name = 'κ', title="Hygroscopicity", ylabel='κ[1]',freq=freq)
    if ('line' in plot):
        line_call(data, slct, x_label="Date", y_label='κ[1]', title='Temporal Variablity in Hygroscopicity')
    input(data)

def Fact_curve(data,freq = 'M', group ='all', ss_choices = ['0.1','0.4','0.7']):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
    group : [str] Time period to generate plots for (default = 'all')
        + 'all' - plot over whole time period given
        + 'year' - generate 1 plot per year if multiple years in data
        + 'season' - generate 1 plot per season
        + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    if len(ss_choices)>1:
        for ss in ss_choices:
            if freq =='S':
                year_vals = np.array([c.split('/')[-1] for c in cols])
                month_vals = np.array([(c.split('_')[-1]).split('/')[0] for c in cols])
                ss_vals = np.array([(c.split('_')[-2]) for c in cols])
                szns= [seasons(float(c)) for c in month_vals]
                select = ['_'.join([s,'/'.join([seasons(float(m)), y])])for s,m, y in zip(ss_vals, month_vals,year_vals)]
                szn_cols = []
                for yr in np.unique(year_vals):
                    for szn in np.unique(szns):
                        s_idx = [i for i,c in enumerate(select) if (szn in c)&(yr in c)&(c.split('_')[0] == ss)]
                        act_col = [cols[i] for i in s_idx]
                        print(act_col)
                        data[f'Fact_{szn}/{yr}'] = data[act_col].mean(axis =1)
                        szn_cols.extend([f'Fact_{szn}/{yr}'])
                plot_data = data[szn_cols]
            if freq =='Y':
                year_vals = np.array([c.split('/')[-1] for c in cols])
                cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
                lowers = [i for i,c in enumerate(cols) if 'lower' in c]
                uppers = [i for i,c in enumerate(cols) if 'upper' in c]
                yr_cols = []
                for yr in np.unique(year_vals):
                    y_idx = [i for i,c in enumerate(year_vals) if (yr in c)&(c.split('_')[0] == ss)]
                    act_col = [cols[i] for i in y_idx]
                    data[f'Fact_{yr}'] = data[act_col].mean(axis =1)
                    yr_cols.extend([f'Fact_{yr}'])
                plot_data = data[yr_cols]
            bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
            slct = [c for c in data.columns.to_numpy() if ('Fact' in c) & (not any(y in c for y in bad_szns))]# & ('autumn' not in c)] 
            plot_data = plot_data[slct]
            plot_data[plot_data>2.1] = np.nan
            plot_data = plot_data.dropna(how='any')
            scat_call(plot_data, slct, yerr_inc= False, name= 'Fact', title = 'Activation Faction',ylabel=r'F$_{act}$[1]',append=ss)
    input(data)

def Fract_curve(data,freq = 'M', group ='all', Kappa =0.1, ss_choices = ['0.1','0.15','0.15','0.4','0.7']):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Activation fraction dataframe
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
    group : [str] Time period to generate plots for (default = 'all')
        + 'all' - plot over whole time period given
        + 'year' - generate 1 plot per year if multiple years in data
        + 'season' - generate 1 plot per season
        + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    if freq =='S':
        szn = ['_'.join([seasons(float(c.split('/')[0])),c.split('/')[-1]]) for c in data.index.to_numpy()]
        data['seasons'] = szn
        print(data)
        data = data.groupby(['seasons']).mean()
        dataf= data.copy()
        ss_vals= sorted({float(re.search(r"\((.*?)\)", col).group(1)) for col in dataf.columns})
        # Build long-form dataframe
        rows = []
        for season in dataf.index:
            for ss in ss_vals:
                rows.append({
                    "season": season,
                    "ss%": ss,
                    "Fact": dataf.loc[season, f"Fact({ss})"],
                    "F>Dcrit": dataf.loc[season, f"F>Dcrit({ss})"]})
        long_df = pd.DataFrame(rows)
        long_df.index
    if freq =='Y':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        yr_cols = []
        for yr in np.unique(year_vals):
            y_idx = [i for i,c in enumerate(year_vals) if (yr in c)&(c.split('_')[0] == ss)]
            act_col = [cols[i] for i in y_idx]
            data[f'Fact_{yr}'] = data[act_col].mean(axis =1)
            yr_cols.extend([f'Fact_{yr}'])
        plot_data = data[yr_cols]
    bad_szns = ['spring_24','autumn_25','summer_26', 'autumn_26']
    plot_data = long_df[~long_df["season"].isin(bad_szns)]
    plot_data = plot_data.dropna(how='any')
    print(np.unique(long_df.season.to_numpy()))
    scat_call_frac_simple(plot_data)
    scat_call_frac_kappa(plot_data)
    scat_call_frac(plot_data,Kappa = Kappa, append=ss)
    for szn_slct in np.unique(long_df.season.to_numpy()):
        szn_data = long_df[long_df["season"]== szn_slct]
        scat_call_frac(szn_data,Kappa = Kappa, append=ss)
    input(data)

'''======Interactive Plotting======+++
Functions used for generating plots:
line_call : calls in the line plot function
line_plot : generates a line plot
scat_call : calls in the scatter plot function
scat_plot : generates a scatter plot
Similar scat plot functionality for the frac curves
with more specificity
+++======Interactive Plotting======'''
def line_call(data, slct, y_label, x_label, title, freq = 'S', ss_vals = ['0.1','0.4','0.7']):
    '''
    Takes in a dataframe of and generates an interactive line plot based on 
    the selected columns.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    none 
    '''
    filt = data[slct] 
    y = []
    dates = []
    leg = []
    for col in slct:
        date = col.split('cutoff_')[-1]
        dates.append(date)
    for ss in ss_vals:
        y.append(filt.filter(items = [ss],axis=0).to_numpy()[0])
        leg.append(f'κ(ss={ss})')
    dates =pd.to_datetime(dates,format='%-d/%-m/%y').to_numpy()
    line_plot(dates,y,leg, y_label=y_label, x_label=x_label, title=title)

def line_plot(x,y,legs, x_label, y_label, title):
    plt.ion()
    fig, ax = plt.subplots()
    fig.autofmt_xdate()
    ax.fmt_xdata = mdates.DateFormatter('%d/%m/%y')
    lines = []
    for i in range(len(y)):
        L, = ax.plot(x,y[i], label = legs[i])
        lines.append(L)
    leg = ax.legend()
    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.set_title(title)
    input('Press enter to exit plot...')
    plt.ioff()

def scat_call(data,slct,freq = 'S',name = r'D$_{crit}$',yerr_inc = True, title = 'Critical diameter', ylabel = 'Diameter [nm]',append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    freq : [str] Frequency of plotting (default = 'S' = season)
    name : [str] variable name for legend (default = "Dcrit")
    yerr_inc : [bool] whether a variable has given error in the Y axis (default = True)
    title : [str] title for plot (default = 'Critical Diameter')
    ylabel : [str] label for y axis (default = 'Diameter [nm]')
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    yerr = []
    fits = []
    clr_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"},}
    colors_month = {'25':{'1':"#701DBE",'2':"#3314BB",'3':"#029B89",
                      '4':"#029D4D",'5':"#00821E",'6':"#A6A902",
                      '7':"#A99002",'8':"#C38500",'9':"#A94502",
                      '10':"#A93A02",'11':"#A90202",'12':"#A90298",}, 
                    '26':{'1':"#9B36F9",'2':"#6D4CFE",'3':"#2DDDC9",
                      '4':"#2ED981",'5':"#2ABF4D",'6':"#DDE028",
                      '7':"#DDC01E",'8':"#FFB005",'9':"#FF6600",
                      '10':"#FE0800",'11':"#FA1937",'12':"#E727D4",}}
    for col in list(slct):
        date = col.split('cutoff_')[-1]
        if 'ideal' in date: 
            color = "#0E1112"
        else:
            if freq == 'M':
                month, yr = date.split('/')
                color = colors_month[yr][month]
            elif (freq == 'd')|(freq == 'w'):
                print(date)
                day,month,yr = date.split('/')
                color = colors_month[yr][month]
            else:
                szn, yr = date.split('/')
                color = colors[szn][yr]
        clr_out.append(color)
        if yerr_inc:
            up = col.replace('cutoff','upper')
            low = col.replace('cutoff','lower')
            upper = abs(data[up].to_numpy()-data[col].to_numpy())
            lower = abs(data[low].to_numpy()-data[col].to_numpy())
            asym = [lower,upper]
            yerr.append(asym)
        else:
            asym = [0.0,0.0]
            yerr.append(asym)
        label = col.split('cuttoff_')[-1]
        y.append(data[col].to_numpy())
        x.append(data.index.to_numpy().astype(float))
        coeff = np.polyfit(data.index.to_numpy().astype(float), data[col].to_numpy(), 5)
        p = np.poly1d(coeff)
        fit = (p(data.index.to_numpy().astype(float)))
        fits.append(fit)
        leg.append(f'{name} {label.split('cutoff_')[-1]}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
    scat_plot(x,y,yerr, fits, leg,clr_out, append=append, title=title, ylabel=ylabel)
    
def scat_plot(x,y, err,fits,legs,colors,append = 0, title = 'Critical diameter', ylabel = 'Diameter [nm]'):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    points = []
    for i in range(len(y)):
        if all([e == 0 for e in err]): 
            ps,caps,bar = ax.errorbar(x[i],y[i], yerr= err[i],ls = '', marker = '.', capsize=0,color= colors[i])
        else:
            ps,caps,bar = ax.errorbar(x[i],y[i], yerr= err[i],ls = '', marker = '.', capsize=3,color= colors[i])
        L = ax.plot(x[i], fits[i],label = legs[i], color = ps.get_color())
        lines.append(L)
        pnt = [ps]
        pnt.extend(list(caps))
        pnt.append(list(bar)[0])
        points.append(pnt)
    leg = ax.legend()
    lined = dict()
    for legline, origline,point in zip(leg.get_lines(), lines, points):       
        legline.set_picker(10)  # 10 pts tolerance
        lined[legline] = origline[0], point
    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline, point = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        for p in point:
            p.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('ss [%]')
    if append ==0:
        ax.set_title(f"{title} of CCN")
    else:
        ax.set_title(f"{title} of CCN for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

'''Frac Plot Functions'''
def scat_call_frac(long_df,name = 'Fract',Kappa =0.1, append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    clr_out = []
    ls_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for season, group in long_df.groupby("season"):
        group = group.sort_values("ss%")
        szn, yr = season.split('_')
        color = colors[szn][yr]
        clr_out.append(color)
        y.append(group["Fact"].to_numpy())
        x.append(group["ss%"].to_numpy().astype(float))
        ls_out.append('-')
        leg.append(r'$F_{act}$ '+ f'{szn}/{yr}')
        clr_out.append(color)
        y.append(group["F>Dcrit"].to_numpy())
        x.append(group["ss%"].to_numpy().astype(float))
        ls_out.append(':')
        leg.append(r'F>$D_{crit}$ '+ f'{szn}/{yr}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
    scat_plot_frac(x,y, leg,clr_out,ls_out, Kappa =Kappa, append=append)
    
def scat_call_frac_kappa(long_df,name = 'Fract',Kappa =0.1):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    clr_out = []
    ls_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for season, group in long_df.groupby("season"):
        group = group.sort_values("ss%")
        szn, yr = season.split('_')
        color = colors[szn][yr]
        clr_out.append(color)
        y.append(group["F>Dcrit"].to_numpy())
        x.append(group["ss%"].to_numpy().astype(float))
        ls_out.append(':')
        leg.append(r'F>$D_{crit}$ '+ f'{szn}/{yr}')
    scat_plot_frac(x,y, leg,clr_out,ls_out, title = f'Fraction of Aerosols larger than Dcrit(κ={Kappa})')

def scat_call_frac_simple(long_df,name = 'Fract'):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    clr_out = []
    ls_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for season, group in long_df.groupby("season"):
        group = group.sort_values("ss%")
        szn, yr = season.split('_')
        color = colors[szn][yr]
        clr_out.append(color)
        y.append(group["Fact"].to_numpy())
        x.append(group["ss%"].to_numpy().astype(float))
        ls_out.append('-')
        leg.append(f'Fact {szn}/{yr}')
    scat_plot_frac(x,y, leg,clr_out,ls_out, title = 'Seasonal Activation Fraction')

def scat_plot_frac(x,y, legs,colors,style, Kappa =0.1, title =0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    for i in range(len(y)):
        L = ax.plot(x[i], y[i], label = legs[i], color = colors[i], ls = style[i])
        lines.append(L)
    leg = ax.legend()
    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):       
        legline.set_picker(10)  # 10 pts tolerance
        lined[legline] = origline[0]
        # input(lined[legline])

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel(r'$Fraction$[1]')
    ax.set_xlabel('ss [%]')
    if title ==0:
        ax.set_title(r"Activation Compared to Aerosols Larger than "+r"$D_{crit}$"+f"(κ={Kappa})")
    else:
        ax.set_title(f"{title}")
    input('Press enter to exit plot...')
    plt.ioff()


"""Calls in the functions"""
if __name__ == '__main__':
    bad_dates = [[pd.to_datetime('10/01/2025 00:00:00'),pd.to_datetime('12/15/2025 00:00:00')],[pd.to_datetime('07/01/2025 00:00:00'),pd.to_datetime('07/24/2025 00:00:00')],[pd.to_datetime('08/06/2025 00:00:00'),pd.to_datetime('08/14/2025 00:00:00')]]
    smps =[r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\2026_SMPS_NumberSizeDist_1hr.csv"]  #list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2026_1hr.csv"]  #r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv"  #[r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv"]#  #list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    master =  r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"
    f = 'W'
    kappa = 0.1
    kappa1 = 0.1
    kappa2 = 0.2
    dataout = r"C:\Users\bensy\Documents\Research\CCN_activation_diameter_test.csv"
    data,ss_cols,diam_cols = comb_files(smps,ccn, freq=f)
    smps_data = data[diam_cols]
    data.to_csv(r"C:\Users\bensy\Documents\Research\check.csv")
    mask = pd.Series(False, index=data.index)
    for date in bad_dates:
        mask |= (data.index >= date[0]) & (data.index <= date[-1])
    data = data[~mask]
    codf = pd.DataFrame()
    kapdf = pd.DataFrame()
    actdf = pd.DataFrame()
    grtr_dict ={'Date':[]}
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        grtr_dict[f'F>Dcrit({ss})'] =[]
        grtr_dict[f'Fact({ss})'] =[]
    for row in range(len(data.index.to_numpy())):
        d = data.iloc[row]
        print(len(d[d.isna()].to_numpy()))
        if len(d[d.isna()].to_numpy())<1:
            cut_offs,kappas = find_cutoff(d,diam_cols,ss_cols)
            print(cut_offs)
            print(kappas)
            act_perc =find_activation(d,diam_cols,ss_cols)
            grtr_dict = find_perc_above(d,diam_cols,ss_cols, grtr_dict, Kappa=kappa)
            Dcrit_ideal1 = {'0.1':[critical_diameter(0.1,kappa=kappa1),critical_diameter(0.1,kappa=kappa1),critical_diameter(0.1,kappa=kappa1)],
                            '0.15':[critical_diameter(0.15,kappa=kappa1),critical_diameter(0.15,kappa=kappa1),critical_diameter(0.15,kappa=kappa1)],
                            '0.25':[critical_diameter(0.25,kappa=kappa1),critical_diameter(0.25,kappa=kappa1),critical_diameter(0.25,kappa=kappa1)],
                            '0.4':[critical_diameter(0.4,kappa=kappa1),critical_diameter(0.4,kappa=kappa1),critical_diameter(0.4,kappa=kappa1)],
                            '0.7':[critical_diameter(0.7,kappa=kappa1),critical_diameter(0.7,kappa=kappa1),critical_diameter(0.7,kappa=kappa1)]}
            Dcrit_ideal2 = {'0.1':[critical_diameter(0.1,kappa=kappa2),critical_diameter(0.1,kappa=kappa2),critical_diameter(0.1,kappa=kappa2)],
                            '0.15':[critical_diameter(0.15,kappa=kappa2),critical_diameter(0.15,kappa=kappa2),critical_diameter(0.15,kappa=kappa2)],
                            '0.25':[critical_diameter(0.25,kappa=kappa2),critical_diameter(0.25,kappa=kappa2),critical_diameter(0.25,kappa=kappa2)],
                            '0.4':[critical_diameter(0.4,kappa=kappa2),critical_diameter(0.4,kappa=kappa2),critical_diameter(0.4,kappa=kappa2)],
                            '0.7':[critical_diameter(0.7,kappa=kappa2),critical_diameter(0.7,kappa=kappa2),critical_diameter(0.7,kappa=kappa2)]}
            if row == 0:
                dt = pd.to_datetime(data.index.to_numpy()[row])
                year = str(dt.year).replace('20','')
                month = dt.month
                day = dt.day
                codf = pd.DataFrame(cut_offs).T
                kapdf= pd.DataFrame(kappas).T
                actdf = pd.DataFrame(act_perc)
                actdf = actdf.set_index('Diameter')
                if f == 'ME':
                    codf.columns = [f'D_cutoff_{month}/{year}',f'D_lower_{month}/{year}',f'D_upper_{month}/{year}']
                    kapdf.columns = [f'Kappa_cutoff_{month}/{year}',f'Kappa_lower_{month}/{year}',f'Kappa_upper_{month}/{year}']
                    actdf = actdf.add_suffix(f'_{month}/{year}')
                else: 
                    codf.columns = [f'D_cutoff_{day}/{month}/{year}',f'D_lower_{day}/{month}/{year}',f'D_upper_{day}/{month}/{year}']
                    kapdf.columns = [f'Kappa_cutoff_{day}/{month}/{year}',f'Kappa_lower_{day}/{month}/{year}',f'Kappa_upper_{day}/{month}/{year}']
                    actdf = actdf.add_suffix(f'_{day}/{month}/{year}')
            else: 
                dt = pd.to_datetime(data.index.to_numpy()[row])
                year = str(dt.year).replace('20','')
                month = dt.month
                day = dt.day
                df = pd.DataFrame(cut_offs).T
                kdf = pd.DataFrame(kappas).T
                adf = pd.DataFrame(act_perc)
                adf = adf.set_index('Diameter')
                if f == 'ME':
                    df.columns = [f'D_cutoff_{month}/{year}',f'D_lower_{month}/{year}',f'D_upper_{month}/{year}']
                    kdf.columns = [f'Kappa_cutoff_{month}/{year}',f'Kappa_lower_{month}/{year}',f'Kappa_upper_{month}/{year}']
                    adf = adf.add_suffix(f'_{month}/{year}')
                else: 
                    df.columns = [f'D_cutoff_{day}/{month}/{year}',f'D_lower_{day}/{month}/{year}',f'D_upper_{day}/{month}/{year}']
                    kdf.columns = [f'Kappa_cutoff_{day}/{month}/{year}',f'Kappa_lower_{day}/{month}/{year}',f'Kappa_upper_{day}/{month}/{year}']
                    adf = adf.add_suffix(f'_{day}/{month}/{year}')
                if (kdf.loc['0.7'][f'Kappa_cutoff_{day}/{month}/{year}']>=0.4).any():
                    print(kdf)
                else:
                    kapdf = pd.merge(kapdf,kdf, left_index=True, right_index=True)
                    codf = pd.merge(codf,df, left_index=True, right_index=True)
                    actdf = pd.merge(actdf, adf, left_index=True, right_index=True)
    # iddf = pd.DataFrame(Dcrit_ideal1).T
    # iddf2 = pd.DataFrame(Dcrit_ideal2).T
    # ideal = [f'D_cutoff_ideal(κ={kappa1})',f'D_lower_ideal(κ={kappa1})',f'D_upper_ideal(κ={kappa1})',f'D_cutoff_ideal(κ={kappa2})',f'D_lower_ideal(κ={kappa2})',f'D_upper_ideal(κ={kappa2})']
    # iddf.columns = [f'D_cutoff_ideal(κ={kappa1})',f'D_lower_ideal(κ={kappa1})',f'D_upper_ideal(κ={kappa1})']
    # iddf2.columns = [f'D_cutoff_ideal(κ={kappa2})',f'D_lower_ideal(κ={kappa2})',f'D_upper_ideal(κ={kappa2})']
    # codf = pd.merge(codf,iddf, left_index=True, right_index=True)
    # codf = pd.merge(codf,iddf2, left_index=True, right_index=True)
    grtr_df = pd.DataFrame(grtr_dict)
    grtr_df = grtr_df.set_index('Date')
    # Fract_curve(grtr_df, freq ='S', Kappa=kappa)
    # Fact_curve(actdf, freq='S')
    # kapdf = kapdf.drop(index='0.1')
    # codf = codf.drop(index='0.1')
    # cut_off_curve(codf, freq ='M')
    kappa_curve(kapdf, freq ='w')
    out = input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        codf.to_csv(out)

