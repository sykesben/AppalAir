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

def line_call(data, slct, append=0):
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
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        y.append(data[ccn_c])
        leg.append(ccn_c)
        if smps_c != 0:
            y.append(data[smps_c])
            leg.append(smps_c)
    line_plot(date,y,leg,append=append)

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
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        y.append(data[ccn_c].to_numpy())
        leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')})')
        y.append(data[smps_c].to_numpy())
        leg.append(f'N{smps_c}')
    hist_plot(y,leg, append=append)

def scat_call(data,slct,append=0):
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
    fit_x =[]
    fit_y = []
    leg = []
    m_all = []
    b_all = []
    r_all = []
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        y.append(data[ccn_c].to_numpy())
        x.append(data[smps_c].to_numpy())
        mask = ~np.isnan(data[smps_c].to_numpy()) & ~np.isnan(data[ccn_c].to_numpy())
        res = linregress(data[smps_c].to_numpy()[mask],data[ccn_c].to_numpy()[mask])
        m,b,r = res.slope,res.intercept,res.rvalue**2
        fit_y.append(m*data[smps_c].to_numpy()+b)
        fit_x.append(data[smps_c].to_numpy())
        rescor = pearsonr(data[ccn_c].to_numpy()[mask],data[smps_c].to_numpy()[mask])
        cor,p = rescor.statistic, rescor.pvalue
        leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')}) vs N{smps_c} | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
        m_all.append(m)
        b_all.append(b)
        r_all.append(cor)
    scat_plot(x,y,fit_x, fit_y,leg,append=append)
    return m_all,b_all,r_all


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
    ax.set_ylabel('N [#/cm^3]')
    ax.set_xlabel('Date')
    if append == 0:
        ax.set_title(f"Comparison of CCN and SMPS data")
    else:
        ax.set_title(f"Comparison of CCN and SMPS data for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

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
    ax.set_ylabel('CCN [#/cm^3]')
    ax.set_xlabel('SMPS [#/cm^3]')
    if append ==0:
        ax.set_title(f"Comparison of CCN and SMPS data")
    else:
        ax.set_title(f"Comparison of CCN and SMPS data for {append}")
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
    ax.set_ylabel('N [#/cm^3]')
    if append == 0: 
        ax.set_title(f"Comparison of CCN and SMPS data")
    else: 
        ax.set_title(f"Comparison of CCN and SMPS data for {append}")
    input('Press enter to exit plot...')
    plt.ioff()