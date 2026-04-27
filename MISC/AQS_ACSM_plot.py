"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between AQS and ACSM
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as LSfit

def line_call(data, slct, append=0):
    '''
    Takes in a dataframe of ACSM and AQS data and generates an interactive line plot based on 
    the selected columns.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined AQS and ACSM data
    slct : [dict of str] relative AQS and ACSM column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    none 
    '''
    date = data.index.to_numpy()
    y = []
    leg = []
    columns_tot = data.columns.to_numpy()
    for AQS_c in list(slct.keys()):
        ACSM_c = slct[AQS_c]
        for c in [c for c in columns_tot if AQS_c in c]:
            y.append(data[c].to_numpy())
            leg.append(f'{c.replace(' [µg/m^3] AQS', '').replace('\n','')}')
        if ACSM_c != 0:
            y.append(data[ACSM_c])
            leg.append(ACSM_c)
    line_plot(date,y,leg,append=append)

def hist_call(data,slct,append=0):
    '''
    Takes in a dataframe of ACSM and AQS data and generates an interactive hist plot based on 
    the selected columns.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined AQS and ACSM data
    slct : [dict of str] relative AQS and ACSM column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    none 
    '''
    y= []
    leg = []
    for AQS_c in list(slct.keys()):
        ACSM_c = slct[AQS_c]
        y.append(data[AQS_c].to_numpy())
        leg.append(f'{AQS_c})')
        y.append(data[ACSM_c].to_numpy())
        leg.append(f'{ACSM_c}')
    hist_plot(y,leg, append=append)


def scat_call(data,slct,append=0):
    '''
    Takes in a dataframe of ACSM and AQS data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined AQS and ACSM data
    slct : [dict of str] relative AQS and ACSM column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    fit_x =[]
    fit_y = []
    leg = []
    m_all = []
    b_all = []
    r_all = []
    for AQS in list(slct.keys()):
        ACSM = slct[AQS]
        y.append(data[AQS].to_numpy())
        x.append(data[ACSM].to_numpy())
        mask = ~np.isnan(data[ACSM].to_numpy()) & ~np.isnan(data[AQS].to_numpy())
        res = linregress(data[ACSM].to_numpy()[mask],data[AQS].to_numpy()[mask])
        m,b,r = res.slope,res.intercept,res.rvalue**2
        fit_y.append(m*data[ACSM].to_numpy()+b)
        fit_x.append(data[ACSM].to_numpy())
        rescor = pearsonr(data[AQS].to_numpy()[mask],data[ACSM].to_numpy()[mask])
        cor,p = rescor.statistic, rescor.pvalue
        leg.append(f'{AQS} vs {ACSM} | {m:.2f}x + {b:.2f} | R2= {r:.4f}')
        m_all.append(m)
        b_all.append(b)
        r_all.append(cor)
    scat_plot(x,y,fit_x, fit_y,leg,append=append)
    return m_all,b_all,r_all

def box_call(data, slct, append=0):
    '''
    Takes in a dataframe of ACSM and AQS data and generates an interactive box plot based on 
    the selected columns.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined AQS and ACSM data
    slct : [dict of str] relative AQS and ACSM column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    none 
    '''
    y= []
    leg = []
    clr =[]
    cols = []
    columns_tot = data.columns.to_numpy()
    for AQS in list(slct.keys()):
        for c in [c for c in columns_tot if AQS in c]:
            cols.append(c)
        cols.append(slct[AQS])
    data = data.dropna(subset = cols,how = 'all')
    for AQS in list(slct.keys()):
        ACSM = slct[AQS]
        for c in [c for c in columns_tot if AQS in c]:
            y.append(data[c].to_numpy())
            if '/total' in c:
                leg.append(f'{c.replace('/total AQS', '\n')}')
                clr.append('skyblue')
            else:
                leg.append(f'{c.replace(' [ug/m3] AQS', '\n')}')
                clr.append('skyblue')
        y.append(data[ACSM].to_numpy())
        if 'total' in ACSM: 
            leg.append(f'{ACSM.replace('/total ', '')}')
        else: 
            leg.append(f'{ACSM.replace('[ug/m3]', '')}')
        clr.append('mediumpurple')
    box_plot(y,leg, clr, append=append)

def scat_plot(x,y,fit_x,fit_y,legs, append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    for i in range(len(y)):
        L, = ax.plot(x[i],y[i], label = legs[i],ls = '', marker = '*')
        lines.append(L)
    fit_lines = []
    leg = ax.legend()
    for i in range(len(fit_y)):
        fit, = ax.plot(fit_x[i],fit_y[i])
        fit_lines.append(fit)
    lined = dict()
    fitlined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline
    for legline, fitline in zip(leg.get_lines(), fit_lines):
        legline.set_picker(5)  # 5 pts tolerance
        fitlined[legline] = fitline

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        fitline = fitlined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        fitline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel('AQS [µg/m^3]')
    ax.set_xlabel('ACSM [µg/m^3]')
    if append ==0:
        ax.set_title(f"Comparison of AQS and ACSM data")
    else:
        ax.set_title(f"Comparison of AQS and ACSM data for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

def box_plot(y,legs, clrs, append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    flier = dict(marker='D', markerfacecolor='orangered', markersize=9,
                  linestyle='none', markeredgecolor='maroon')
    median= dict(color = 'maroon', linewidth = 3)
    bplot = ax.boxplot(y,
                flierprops=flier,
                medianprops=median,
                patch_artist=True, # color plots
                tick_labels=legs) # will be used to label x-ticks)  
    # fill with colors
    for patch, color in zip(bplot['boxes'], clrs):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
    if append == 0: 
        ax.set_title(f"Comparison of AQS and ACSM data")
    else: 
        ax.set_title(f"Comparison of AQS and ACSM data for {append}")
    ax.set_ylabel('X[µg/total]')
    input('Press enter to exit plot...')
    plt.ioff()

def line_plot(x,y,legs, append = 0):
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
    ax.set_ylabel('X [µg/m^3]')
    ax.set_xlabel('Date')
    if append == 0:
        ax.set_title(f"Comparison of AQS and ACSM data")
    else:
        ax.set_title(f"Comparison of AQS and ACSM data for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

def hist_plot(y,legs, append = 0):
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
    ax.set_ylabel('X [µg/m^3]')
    if append == 0: 
        ax.set_title(f"Comparison of AQS and ACSM data")
    else: 
        ax.set_title(f"Comparison of AQS and ACSM data for {append}")
    input('Press enter to exit plot...')
    plt.ioff()
