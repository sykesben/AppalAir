"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as LSfit
plt.rcParams['font.size'] = 20
# plt.rcParams['font.weight'] = 'bold'

def line_call(data, slct, y_label, x_label, title):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive line plot based on 
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
    date = data.index.to_numpy()
    y = []
    leg = []
    for var_1 in list(slct.keys()):
        var_2 = slct[var_1]
        y.append(data[var_1])
        leg.append(var_1)
        if var_2 != 0:
            y.append(data[var_2])
            leg.append(var_2)
    line_plot(date,y,leg, y_label=y_label, x_label=x_label, title=title)

def hist_call(data,slct,append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive hist plot based on 
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
    y= []
    leg = []
    for var_1 in list(slct.keys()):
        var_2 = slct[var_1]
        y.append(data[var_1].to_numpy())
        leg.append(f'{var_1}')
        y.append(data[var_2].to_numpy())
        leg.append(f'{var_2}')
    hist_plot(y,leg, append=append)

def scat_call(data, slct, x_label, y_label, title, split=0, verbose = True, single = True):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])
    split : [dict or int] list of values to split on (default = 0[no split])

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
    if split == 0:
        for var_1 in list(slct.keys()):
            var_2 = slct[var_1]
            y.append(data[var_1].to_numpy())
            x.append(data[var_2].to_numpy())
            mask = ~np.isnan(data[var_2].to_numpy()) & ~np.isnan(data[var_1].to_numpy())
            res = linregress(data[var_2].to_numpy()[mask],data[var_1].to_numpy()[mask])
            m,b,r = res.slope,res.intercept,res.rvalue**2
            fit_y.append(m*data[var_2].to_numpy()+b)
            fit_x.append(data[var_2].to_numpy())
            rescor = pearsonr(data[var_1].to_numpy()[mask],data[var_2].to_numpy()[mask])
            cor,p = rescor.statistic, rescor.pvalue
            if single: 
                if verbose: 
                    leg.append(f'{m:.2e}x + {b:.2f} | corr = {cor:.2f}')
                else: 
                    leg.append(f'{m:.2e}x + {b:.2f}')
            if verbose: 
                leg.append(f'{var_1} vs {var_2} | {m:.2e}x + {b:.2f} | corr = {cor:.2f}')
            else: 
                leg.append(f'{var_1} vs {var_2} ({m:.2e}x + {b:.2f})')
            m_all.append(m)
            b_all.append(b)
            r_all.append(cor)
    else:
        key = list(split.keys())[0]
        input(key)
        for val in list(split.values())[0]:
            print(val)
            append = val[0]
            vals = val[1]
            sub_indices = data.index[data[key].isin(vals)] 
            subdata = data.drop(sub_indices)
            for var_1 in list(slct.keys()):
                var_2 = slct[var_1]
                y.append(subdata[var_1].to_numpy())
                x.append(subdata[var_2].to_numpy())
                mask = ~np.isnan(subdata[var_2].to_numpy()) & ~np.isnan(subdata[var_1].to_numpy())
                res = linregress(subdata[var_2].to_numpy()[mask],subdata[var_1].to_numpy()[mask])
                m,b,r = res.slope,res.intercept,res.rvalue**2
                fit_y.append(m*subdata[var_2].to_numpy()+b)
                fit_x.append(subdata[var_2].to_numpy())
                rescor = pearsonr(subdata[var_1].to_numpy()[mask],subdata[var_2].to_numpy()[mask])
                cor,p = rescor.statistic, rescor.pvalue
                leg.append(f'{var_1} vs {var_2} at {append} | corr = {cor:.2f}')
                m_all.append(m)
                b_all.append(b)
                r_all.append(cor)
    scat_plot(x,y,fit_x, fit_y, leg, x_label, y_label, title)
    return m_all,b_all,r_all

def box_call(data, slct, y_label, title):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive box plot based on 
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
    y= []
    leg = []
    clr =[]
    cols = []
    for var_1 in list(slct.keys()):
        cols.append(var_1)
        cols.append(slct[var_1])
    data = data.dropna(subset = cols)
    for var_1 in list(slct.keys()):
        var_2 = slct[var_1]
        y.append(data[var_1].to_numpy())
        leg.append(f'{var_1})')
        clr.append('skyblue')
        y.append(data[var_2].to_numpy())
        leg.append(f'{var_2}')
        clr.append('mediumpurple')
    box_plot(y,leg, clr, y_label=y_label, title=title)

def line_plot(x,y,legs, x_label, y_label, title):
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
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.set_title(title)
    input('Press enter to exit plot...')
    plt.ioff()

def scat_plot(x,y,fit_x,fit_y,legs, x_label, y_label, title):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    for i in range(len(y)):
        L, = ax.plot(x[i],y[i], label = legs[i],ls = '', marker = '*')
        lines.append(L)
    fit_lines = []
    leg = ax.legend()
    for i in range(len(fit_y)):
        fit, = ax.plot(fit_x[i],fit_y[i], color = lines[i].get_color())
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
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.set_title(title)
    input('Press enter to exit plot...')
    plt.ioff()

def hist_plot(y,legs, y_label, title):
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
    ax.set_ylabel(y_label)
    ax.set_title(title)
    input('Press enter to exit plot...')
    plt.ioff()

def box_plot(y,legs, clrs, y_label, title):
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
    ax.set_title(title)
    ax.set_ylabel(y_label)
    input('Press enter to exit plot...')
    plt.ioff()