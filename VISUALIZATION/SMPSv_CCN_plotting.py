"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
from SMPSvCCNplot_gen import line_call, hist_call,scat_call, box_call,chem_line_call,chem_scat_call, cor_scat_call, cor_box_call, cor_line_call
from pull_in import *

def plot_gen(data, ss2nm, mode = 0,vars = ['ss'], date = 0, group ='all', thresh = 0, cormode = False):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    vars : [list of str] columns to use while plotting (default = ['ss'])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)
    cormode : [bool] for comparing between corrections (default = False)

    Returns
    ++++++++++
    none 
    '''
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
    ccn_cols =[]
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
    # if cormode:
    #     ccn_cors = [list(col.split('setpt'))[-1] for col in cols if "N(cm-3)_cor_" in col]
    #     ccn_cols = []
    #     choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
    val =0
    title = 0
    if 'Fact' in vars:
        cpc_col = 'F_act_CPC(ss=0.7%)'
        smps_col ='F_act_SMPS(ss=0.7%)'
        slct[cpc_col] = smps_col
        val = 'Fact[1]'
        title = 'Activation Fraction at 0.7%ss'
    if 'counts' in vars:
        ccn_col = f'N(cm-3)_cor_setpt0.7'
        cpc_col = f'CPC N[cm-3]'
        smps_col = f'Total Concentration (#/cm³)'
        slct[ccn_col] = [smps_col,cpc_col]
        ccn_cols.append(ccn_col)
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        slct[f'Total Concentration (#/cm³)'] = f'CPC N[cm-3]'
        if choice == 'all':
            for ss in ss_vals:
                ccn_col = f'N(cm-3)_cor_setpt{ss}'
                smps_col = f'{ss2nm[ss]}'
                ccn_cols.append(ccn_col)
                slct[ccn_col] = smps_col
        else:
            ccn_col = f'N(cm-3)_cor_setpt{choice}'
            smps_col = f'{ss2nm[choice]}'
            slct[ccn_col] = smps_col
            ccn_cols.append(ccn_col)
    for c in ccn_cols:
        drop_indices = data[data[c]<thresh].index   
        data = data.drop(drop_indices)
    if ('Q' in vars) & (mode == 'line'): slct['Q(lpm)_sample'] = 0
    if ('T' in vars) & (mode == 'line'): slct['T(C)_sample'] = 0 
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            if cormode:
                cor_line_call(data,slct)
            else:
                line_call(data,slct, val=val, title=title)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            if cormode:
                cor_scat_call(data,slct)
            else:
                scat_call(data, slct)
        elif mode == 'hist': #histogram
            hist_call(data,slct)
        elif mode == 'box':# box and whisker
            if cormode:
                cor_box_call(data,slct)
            else: 
                box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                if cormode:
                    cor_line_call(Ydata,slct,append)
                else:
                    line_call(Ydata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                if cormode:
                    m,b,cor = cor_scat_call(Ydata, slct,append)
                else:
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
                if cormode: 
                    cor_box_call(Ydata,slct,append)
                else:
                    box_call(Ydata,slct,append)
    elif group == 'month':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    if cormode:
                        cor_line_call(Mdata,slct,append)
                    else:
                        line_call(Mdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    if cormode:
                        m,b,cor = cor_scat_call(Mdata, slct, append)
                    else:
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
                    if cormode:
                        cor_box_call(Mdata,slct,append)
                    else:
                        box_call(Mdata,slct,append)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    if cormode:
                        cor_line_call(Sdata,slct,append)
                    else:
                        line_call(Sdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    if cormode:
                        m,b,cor = scat_call(Sdata,slct, append)
                    else:
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
                    if cormode:
                        cor_box_call(Sdata,slct,append)
                    else: 
                        box_call(Sdata,slct,append)

def chem_plot_gen(data, mode = 'line',vars = ['ss'], date = 0, group ='all', thresh = 0, chem = 'org'):
    '''
    Takes in a dataframe of SMPS and CCN and chemistry data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = line)
    vars : [list of str] columns to use while plotting (default = ['ss'])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)
    chem : [str] Species to compare data to (default = 'org')

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
    comp_cols =[]
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
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        if choice == 'all':
            for ss in ss_vals:
                ccn_col = f'N(cm-3)_cor_setpt{ss}'
                smps_col = f'>{ss2nm[ss]}nm'
                new_col = f'CCN and SMPS %dev at ss={ss}'
                A = data[ccn_col].to_numpy()
                B = data[smps_col].to_numpy()
                perc_dev = (np.abs(A-B))/((A+B)/2)*100
                data[new_col] = perc_dev
                comp_cols.append(new_col)
        else:
            ccn_col = f'N(cm-3)_cor_setpt{choice}'
            smps_col = f'>{ss2nm[choice]}nm'
            new_col = f'CCN and SMPS %dev at ss={choice}'
            A = data[ccn_col].to_numpy()
            B = data[smps_col].to_numpy()
            perc_dev = (np.abs(A-B))/((A+B)/2)*100
            data[new_col] = perc_dev
            comp_cols.append(new_col)
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            chem_line_call(data,comp_cols, chem=chem)
        elif mode =='scat':
            chem_scat_call(data, comp_cols, chem=chem)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                chem_line_call(Ydata,comp_cols,append=append, chem=chem)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                m,b,cor = chem_scat_call(Ydata,comp_cols,append=append, chem=chem)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    chem_line_call(Sdata,comp_cols,append=append, chem=chem)
        elif mode == 'scat': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    chem_scat_call(Sdata,comp_cols,append=append, chem=chem)

def plot_gen_corr(data, mode = 0,vars = ['ss'], date = 0, group ='all', thresh = 0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    vars : [list of str] columns to use while plotting (default = ['ss'])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)

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
    smps_cols =[]
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
        mode = input('What style plot would you like to generate? (line, scat, hist, box) ')
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        if choice == 'all':
            for ss in ss_vals:
                ccn_cols = [col for col in data.columns.to_numpy() if f'N(cm-3)_cor_setpt{ss}' in col]
                smps_col = f'>{ss2nm[ss]}nm'
                smps_cols.append(smps_col)
                slct[smps_col] = ccn_cols
        else:
            ccn_cols = [col for col in data.columns.to_numpy() if f'N(cm-3)_cor_setpt{choice}' in col]
            smps_col = f'>{ss2nm[choice]}nm'
            smps_cols.append(smps_col)
            slct[smps_col] = ccn_cols
    for c in ccn_cols:
        drop_indices = data[data[c]<thresh].index   
        data = data.drop(drop_indices)
    if ('Q' in vars) & (mode == 'line'): slct['Q(lpm)_sample'] = 0
    if ('T' in vars) & (mode == 'line'): slct['T(C)_sample'] = 0 
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            line_call(data,slct)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            scat_call(data, slct)
        elif mode == 'hist': #histogram
            hist_call(data,slct)
        elif mode == 'box':# box and whisker
            box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                line_call(Ydata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
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
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    line_call(Mdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
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
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    line_call(Sdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
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

if __name__ == '__main__':
    #list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    cpc = r"C:\Users\bensy\Documents\Research\app_CPC.csv"
    smps =[r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\2026_SMPS_NumberSizeDist_1hr.csv"]  #list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2026_1hr.csv"]     #list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    # AQS = [r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_avery_2024.csv",r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_avery_2025.csv",r"C:\Users\bensy\Documents\Research\AQS_Processed\AQS_avery_2025.csv"]
    master =  r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"

    data,ss2nm = comb_files(smps, ccn, cpc,freq='d',kappa =0.1)
    # data_chem = comb_files(smps,ccn,cpc,freq="W")
    input(data)
    out = r"C:\Users\bensy\Documents\Research\CCN_vSMPS.csv"#input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        data.to_csv(out)
    # bad_dates = ['8/15/2025 00:00:00','12/01/2025 00:00:00']
    plot_gen(data,ss2nm,mode = 'line',thresh=5,group='all',vars=['ss'])
    # chem_plot_gen(data_chem,mode = 'line',chem = 'org')