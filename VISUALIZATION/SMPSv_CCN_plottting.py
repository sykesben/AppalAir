"""
Date: 1/28/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
import os
import datetime as dt 
from scipy.stats import linregress
import matplotlib.pyplot as plt
from pathlib import Path
from os.path import expanduser 
from scipy.optimize import least_squares as LSfit


def smps_data(files,freq='d',ss = [0.1,0.7]):
    '''
    Takes in a list of smps files and filters out the particle size concentration depending
    on comparable ss% values.
    ----------

    Paramaters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame
    ss : [list of floats] ss% set points from CCN

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
        if i == 0:
            smps= file
        else:
            smps = pd.concat([smps,file])
    smps.index = pd.to_datetime(smps.index)
    numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]
    smps = smps[numsmps]
    ss2nm = {0.1:200, 0.6:100, 0.7:80}
    cols = []
    for s in ss:
        col = f'>{ss2nm[s]}nm'
        cols.append(col)
        gr = [n for n in numsmps if float(n) >ss2nm[s]]
        smps[col] = smps[gr].mean(axis=1)  
    smps = smps.resample(freq).mean()
    smps.index.names = ['Date']
    return smps,cols

def ccn_data(files, freq ='d', ss = [0.1,0.7]):
    '''
    Takes in a list of CCN files and returns a processed dataframe with important columns 
    for plotting or further analysis
    ----------

    Paramaters
    ++++++++++
    files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrame (default = 'd')
    ss : [list of floats] ss% set points from CCN (default = [0.1,0.7])

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
            file=file.set_index('Datetime UTC') #Set index
        if i == 0:
            ccn = file
        else:
            ccn = pd.concat([ccn,file])
    ccn.index = pd.to_datetime(ccn.index)
    cols = ['T(C)_inlet','T1(C)','T(C)_sample','T(C)_OPC','T(C)_nafion','Q(lpm)_sample','Q(lpm)_sheath','P(hPA)_sample']
    for s in ss:
        cols.append(f'N(cm-3)_avg_setpt{s}')
        cols.append(f'N(cm-3)_cor_setpt{s}')
        cols.append(f'ss(%)_calc_setpt{s}')
    ccn = ccn[cols]
    ccn = ccn.resample(freq).mean()
    ccn.index.names = ['Date']
    return ccn,cols

def comb_files(smps_files,ccn_files, ss = [0.1,0.7], freq = 'd'):
    '''
    Takes in a list of CCN and SMPS files and returns a combined dataframe with important columns 
    from both for plotting or further analysis
    ----------

    Paramaters
    ++++++++++
    smps_files : [list of str] Paths to SMPS files
    ccn_files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrames (default = 'd')
    ss : [list of floats] ss% set points from CCN (default = [0.1,0.7])

    Returns
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data from all inputted files
    '''
    ccn, ccn_cols = ccn_data(ccn_files,freq,ss)
    smps, smps_cols = smps_data(smps_files,freq,ss)
    data = pd.merge(ccn[ccn_cols],smps[smps_cols],left_index = True, right_index = True)
    return data

def plot_gen(data, mode = 0,vars = ['ss'], date = 0, group ='all'):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Paramaters
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
    input(data.season)
    slct = {}
    cols = data.columns.to_numpy()
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
                slct[ccn_col] = smps_col
        else:
            ccn_col = f'N(cm-3)_cor_setpt{choice}'
            smps_col = f'>{ss2nm[choice]}nm'
            slct[ccn_col] = smps_col
    if ('Q' in vars) & (mode == 'line'): slct['Q(lpm)_sample'] = 0
    if ('T' in vars) & (mode == 'line'): slct['T(C)_sample'] = 0 
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            date = data.index.to_numpy()
            y = []
            leg = []
            for ccn_c in list(slct.keys()):
                smps_c = slct[ccn_c]
                y.append(data[ccn_c])
                leg.append(ccn_c)
                if smps_c != 0:
                    y.append(data[smps_c])
                    leg.append(smps_c)
            line_plot(date,y,leg)

        elif mode == 'scat': # scatter plot SMPS vs CCN
            x = []
            y = []
            leg = []
            for ccn_c in list(slct.keys()):
                smps_c = slct[ccn_c]
                y.append(data[ccn_c].to_numpy())
                x.append(data[smps_c].to_numpy())
                mask = ~np.isnan(data[smps_c].to_numpy()) & ~np.isnan(data[ccn_c].to_numpy())
                res = linregress(data[smps_c].to_numpy()[mask],data[ccn_c].to_numpy()[mask])
                m,b,r = res.slope,res.intercept,res.rvalue**2
                leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')}) vs N{smps_c} | {m:.2f}x + {b:.2f} | R2= {r:.4f} ')
            scat_plot(x,y,leg)

        elif mode == 'hist': #histogram
            y= []
            leg = []
            for ccn_c in list(slct.keys()):
                smps_c = slct[ccn_c]
                y.append(data[ccn_c].to_numpy())
                leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')})')
                y.append(data[smps_c].to_numpy())
                leg.append(f'N{smps_c}')
            hist_plot(y,leg)
    if group == 'year':
        years ={}
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                date = Ydata.index.to_numpy()
                y = []
                leg = []
                for ccn_c in list(slct.keys()):
                    smps_c = slct[ccn_c]
                    y.append(Ydata[ccn_c])
                    leg.append(ccn_c)
                    if smps_c != 0:
                        y.append(Ydata[smps_c])
                        leg.append(smps_c)
                line_plot(date,y,leg)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                x = []
                y = []
                leg = []
                for ccn_c in list(slct.keys()):
                    smps_c = slct[ccn_c]
                    y.append(Ydata[ccn_c].to_numpy())
                    x.append(Ydata[smps_c].to_numpy())
                    mask = ~np.isnan(Ydata[smps_c].to_numpy()) & ~np.isnan(Ydata[ccn_c].to_numpy())
                    res = linregress(Ydata[smps_c].to_numpy()[mask],Ydata[ccn_c].to_numpy()[mask])
                    m,b,r = res.slope,res.intercept,res.rvalue**2
                    leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')}) vs N{smps_c} | {m:.2f}x + {b:.2f} | R2= {r:.4f} ')
                scat_plot(x,y,leg)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                y= []
                leg = []
                for ccn_c in list(slct.keys()):
                    smps_c = slct[ccn_c]
                    y.append(Ydata[ccn_c].to_numpy())
                    leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')})')
                    y.append(Ydata[smps_c].to_numpy())
                    leg.append(f'N{smps_c}')
                hist_plot(y,leg)
            


def line_plot(x,y,legs):
    plt.ion()
    fig, ax = plt.subplots()
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
    ax.set_ylabel('N [#/cm^3]')
    ax.set_xlabel('Date')
    ax.set_title(f"Comparison of CCN and SMPS data")
    input('Press enter to exit plot...')
    plt.ioff()

def scat_plot(x,y,legs):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    for i in range(len(y)):
        L, = ax.plot(x[i],y[i], label = legs[i],ls = '', marker = '*')
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
    ax.set_ylabel('CCN [#/cm^3]')
    ax.set_xlabel('SMPS [#/cm^3]')
    ax.set_title(f"Comparison of CCN and SMPS data")
    input('Press enter to exit plot...')
    plt.ioff()

def hist_plot(y,legs):
    plt.ion()
    fig, ax = plt.subplots()
    art = []
    for i in range(len(y)):
        # input(ax.hist(y[i], label = legs[i]),alpha=0.5)
        n,b,a = ax.hist(y[i], label = legs[i],alpha=0.5)
        art.append(a)
    leg = ax.legend()
    lined = dict()
    for legline, origline in zip(leg.get_lines(), art):
        input(legline)
        legline.set_picker(5)  # 5 pts tolerance
        input(origline)
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
    ax.set_ylabel('N [#/cm^3]')
    ax.set_title(f"Comparison of CCN and SMPS data")
    input('Press enter to exit plot...')
    plt.ioff()


if __name__ == '__main__':
    smps = [r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv"]#list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv"]#list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    data = comb_files(smps,ccn)
    out = input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        data.to_csv(out)
    plot_gen(data, group ='Seasonal')