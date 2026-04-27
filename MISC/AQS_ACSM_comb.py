"""
Date: 3/18/26
Author: Ben Sykes
Purpose: generate plots between AQS and ACSM datasets
"""

import pandas as pd
import numpy as np
from os.path import basename, join, dirname,expanduser
pd.set_option('mode.chained_assignment', None)
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as LSfit
from AQS_ACSM_plot import line_call, hist_call,scat_call, box_call


def AQS_CSVs_for_Reindexing(files,freq='W'):
    aqs = pd.DataFrame()
    for i in range(len(files)): #read in AQS files and combine
        f = files[i]
        file =pd.read_csv(f) #read in AQS file
        file=file.set_index("Date(UTC)") #Set index
        if i == 0:
            aqs = file
        else:
            aqs = pd.concat([aqs,file])
    aqs.index = pd.to_datetime(aqs.index)
    specs = ['NH4/total','EC/total','OC/total','SO4/total','NO3/total', 'Org/total']
    locs = aqs['Location'].unique()
    aqs['Org/total'] = aqs['OC/total']*2
    AQS_tot = pd.DataFrame()
    for loc in locs:
        aqs_locs = aqs[aqs['Location']== loc]
        aqs_locs = aqs_locs[specs]
        state_dict = {'North Carolina': 'NC',
                      'Georgia': 'GA',
                      'Kentucky': 'KY',
                      'South Carolina': 'SC',
                      'Tennessee': 'TN',
                      'Virginia':'VA',
                      'West Virginia': 'WV'}
        loc = loc.replace('county', '')
        name = list(loc.split(','))
        name[-1] = state_dict[name[-1]]
        loc = '\n'.join(name)
        aqs_locs = aqs_locs.add_suffix(f' {loc}')
        if loc == locs[0]:
            AQS_tot = aqs_locs
        else: 
            AQS_tot = pd.concat([AQS_tot,aqs_locs],axis=1)
    AQS_tot.columns = AQS_tot.columns.str.replace('/total', '/total AQS')
    AQS_tot = AQS_tot.resample(freq).apply(np.nanmean)
    AQS_tot = AQS_tot.dropna(how = 'all')
    return AQS_tot,specs

def PM25_data(files,freq='W'):
    aqs = pd.DataFrame()
    for i in range(len(files)): #read in AQS files and combine
        f = files[i]
        file =pd.read_csv(f) #read in AQS file
        file=file.set_index("Date(UTC)") #Set index
        if i == 0:
            aqs = file
        else:
            aqs = pd.concat([aqs,file])
    aqs.index = pd.to_datetime(aqs.index)
    specs = ['PM2.5 [ug/m3 STP]']
    aqs = aqs[specs]
    aqs.columns = aqs.columns.str.replace('[ug/m3 STP]', '[ug/m3] AQS')
    aqs = aqs.resample(freq).apply(np.nanmean)
    aqs = aqs.dropna()
    return aqs,specs

def master_data(f,freq='W'):
    master=pd.read_csv(f) #read in AQS file
    master=master.set_index("Local time (UTC-5)") #Set index
    # input(master.columns)
    master['Date(UTC)'] = pd.to_datetime(master.index) + pd.Timedelta(hours=5)
    specs = ['NH4_11000','SO4_11000','NO3_11000','Org_11000','SO4/total','org/total','1hrMC_µg/m3']
    master=master.set_index("Date(UTC)")
    master = master[specs]
    master.columns = master.columns.str.replace('_11000', ' [ug/m3] ACSM')
    master.columns = master.columns.str.replace('/total', '/total ACSM')
    master.columns = master.columns.str.replace('org', 'Org')
    master.columns = master.columns.str.replace('1hrMC_µg/m3', 'PM2.5 [ug/m3] ACSM')
    master = master.resample(freq).apply(np.nanmean)
    master = master.dropna()
    return master,specs

def plot_gen(data, mode = 0,vars = ['spec'], date = 0, group ='all'):
    '''
    Takes in a dataframe of AQS and ACSM data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined AQS and ACSM data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    vars : [list of str] columns to use while plotting (default = [spec], or [pm25])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    ss2nm = {'0.1':200, '0.6':100, '0.7':80}
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    data['year'] = data.index.year.to_numpy()
    data['month']= data.index.month.to_numpy()
    data['season'] = [seasons(n) for n in data.index.month.to_numpy()]
    slct = {}
    cols = data.columns.to_numpy()
    AQS_cols =[]
    if date != 0: #if date is passed, split the data using the passed date range
        if isinstance(date, list):
            if len(date) > 1:
                date0 = pd.to_datetime(date[0])
                date1 = pd.to_datetime(date[1])
                data = data.loc[date0:date1]
            else:
                date0 = pd.to_datetime(date[0])
                data = data.loc[date0:]
        elif isinstance(date, str):
            date0 = pd.to_datetime(date)
            data = data.loc[date0:]
    if mode == 0: #if default value used and no mode passed, user input mode
        mode = input('What style plot would you like to generate? (line, scat, hist)')
    if 'spec' in vars: #if speciated values used for plotting
        chem_vals = ['NH4' ,'SO4' ,'NO3', 'Org']
        choice = input(f'Which speciated data value would you like to see? ({', '.join(chem_vals)}, or all) ')
        if choice == 'all':
            for chem in chem_vals:
                AQS_col = f'{chem}/total AQS'
                ACSM_col = f'{chem}/total ACSM'
                slct[AQS_col] = ACSM_col
                AQS_cols.append(AQS_col)
        else:
            AQS_col = f'{choice}/total AQS'
            ACSM_col = f'{choice}/total ACSM'
            slct[AQS_col] = ACSM_col
            AQS_cols.append(AQS_col)
    elif 'pm25' in vars: #if speciated values used for plotting
        AQS_col = f'PM2.5 [ug/m3] AQS'
        ACSM_col = f'PM2.5 [ug/m3] ACSM'
        AQS_cols.append(AQS_col)
        slct[AQS_col] = ACSM_col
    if group =='all':
        if mode == 'line': ## Line plot ACSM and AQS vs date
            line_call(data,slct)
        elif mode == 'scat': # scatter plot ACSM vs AQS
            scat_call(data, slct)
        elif mode == 'hist': #histogram
            hist_call(data,slct)
        elif mode == 'box':# box and whisker
            box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot ACSM and AQS vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                line_call(Ydata,slct,append)
        elif mode == 'scat': # scatter plot ACSM vs AQS
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                m,b,cor = scat_call(Ydata, slct,append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                hist_call(Ydata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                box_call(Ydata,slct,append)
    elif group == 'month':
        if mode == 'line': ## Line plot ACSM and AQS vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    line_call(Mdata,slct,append)
        elif mode == 'scat': # scatter plot ACSM vs AQS
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    m,b,cor = scat_call(Mdata, slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    hist_call(Mdata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    box_call(Mdata,slct,append)
    elif group == 'season':
        if mode == 'line': ## Line plot ACSM and AQS vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    line_call(Sdata,slct,append)
        elif mode == 'scat': # scatter plot ACSM vs AQS
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    m,b,cor = scat_call(Sdata,slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    hist_call(Sdata,slct,append)
        elif mode == 'box': #box and whisker plots
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    box_call(Sdata,slct,append)

f = [r"C:\Users\bensy\Documents\Research\AQS_CSVs_for_Reindexing\AQS_combined_speciated_2024.csv",r"C:\Users\bensy\Documents\Research\AQS_CSVs_for_Reindexing\AQS_combined_Speciated_2025.csv"]
aqs,specs = AQS_CSVs_for_Reindexing(f)
f2 = [r"C:\Users\bensy\Documents\Research\AQS_CSVs_for_Reindexing\AQS_combined_PM25_2024.csv",r"C:\Users\bensy\Documents\Research\AQS_CSVs_for_Reindexing\AQS_combined_PM25_2025.csv"]
pm25,specs = PM25_data(f2)
f3 = r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"
master, specs = master_data(f3)
data = pd.merge(aqs,master,left_index = True, right_index = True)
data = pd.merge(data,pm25,left_index = True, right_index = True)
data = data.dropna(axis =1, how = 'all')
data = data.ffill().bfill()
# data = data.fillna(0) 
# input(data)
data.to_csv(expanduser("~/Documents/Research/Chemistry_comparison_daily.csv"))
plot_gen(data,mode ='box',vars=['spec'])
