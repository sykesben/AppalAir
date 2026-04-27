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

def chem_line_call(data, comp, chem ='org', append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive line plot based on 
    the selected columns.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    comp : [list of str] comparison CCN and SMPS column names for processing
    chem : [str] chemical species to use during comparison 
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    none 
    '''
    date = data.index.to_numpy()
    y = []
    leg = []
    chem_trend = f'{chem}/total'
    y.append(data[chem_trend].to_numpy()*100)
    leg.append(f'{chem} %')
    for c in comp:
        y.append(data[c])
        leg.append(c)
    chem_line_plot(date,y,leg,append=append)

def chem_scat_call(data,comp,chem ='org',append=0):
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
    chem_trend = f'{chem}/total'
    for dev in comp:
        y.append(data[dev].to_numpy())
        x.append(data[chem_trend].to_numpy()*100)
        mask = ~np.isnan(data[chem_trend].to_numpy()) & ~np.isnan(data[dev].to_numpy())
        res = linregress(data[chem_trend].to_numpy()[mask],data[dev].to_numpy()[mask])
        m,b,r = res.slope,res.intercept,res.rvalue**2
        fit_y.append(m*data[chem_trend].to_numpy()+b)
        fit_x.append(data[chem_trend].to_numpy())
        rescor = pearsonr(data[dev].to_numpy()[mask],data[chem_trend].to_numpy()[mask])
        cor,p = rescor.statistic, rescor.pvalue
        leg.append(f'{dev} vs {chem}% | corr = {cor:.2f}')
        m_all.append(m)
        b_all.append(b)
        r_all.append(cor)
    chem_scat_plot(x,y,fit_x, fit_y,leg,append=append)
    return m_all,b_all,r_all

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

def cor_line_call(data, slct, append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive line plot based on 
    the selected columns. Comparing different corr mode to base data.
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
    columns = data.columns.to_numpy()
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        for col in [c for c in columns if (f'{ccn_c}_' in c) |(c ==ccn_c)]:
            y.append(data[col].to_numpy())
            leg.append(col)
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
        fit_y.append(1*data[smps_c].to_numpy())
        fit_x.append(data[smps_c].to_numpy())
        rescor = pearsonr(data[ccn_c].to_numpy()[mask],data[smps_c].to_numpy()[mask])
        cor,p = rescor.statistic, rescor.pvalue
        leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')}) vs N{smps_c}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
        m_all.append(m)
        b_all.append(b)
        r_all.append(cor)
    scat_plot(x,y,fit_x, fit_y,leg,append=append)
    return m_all,b_all,r_all

def cor_scat_call(data,slct,append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode. Comparing different corr mode
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
    cols= data.columns.to_numpy()
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        x.append(data[smps_c].to_numpy())
        for col in [c for c in cols if ccn_c in c]:
            y.append(data[ccn_c].to_numpy())
        
        mask = ~np.isnan(data[smps_c].to_numpy()) & ~np.isnan(data[ccn_c].to_numpy())
        res = linregress(data[smps_c].to_numpy()[mask],data[ccn_c].to_numpy()[mask])
        m,b,r = res.slope,res.intercept,res.rvalue**2
        fit_y.append(1*data[smps_c].to_numpy())
        fit_x.append(data[smps_c].to_numpy())
        rescor = pearsonr(data[ccn_c].to_numpy()[mask],data[smps_c].to_numpy()[mask])
        cor,p = rescor.statistic, rescor.pvalue
        leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')}) vs N{smps_c}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
        m_all.append(m)
        b_all.append(b)
        r_all.append(cor)
    scat_plot(x,y,fit_x, fit_y,leg,append=append)
    return m_all,b_all,r_all

def box_call_corr(data, slct, append=0):
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
    for ccn_c in list(slct.keys()):
        cols.append(ccn_c)
        cols.append(slct[ccn_c])
    data = data.dropna(subset = cols)
    for smps_c in list(slct.keys()):
        ccn_cs = slct[smps_c]
        for ccn_c in ccn_cs:
            y.append(data[ccn_c].to_numpy())
            leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')})')
            clr.append('skyblue')
        y.append(data[smps_c].to_numpy())
        leg.append(f'N{smps_c}')
        clr.append('mediumpurple')
    box_plot(y,leg, clr, append=append)

def box_call(data, slct, append=0):
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
    for ccn_c in list(slct.keys()):
        cols.append(ccn_c)
        cols.append(slct[ccn_c])
    data = data.dropna(subset = cols)
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        y.append(data[ccn_c].to_numpy())
        leg.append(f'{ccn_c.replace('cm-3)_cor_setpt', 'ss%=')})')
        clr.append('skyblue')
        y.append(data[smps_c].to_numpy())
        leg.append(f'N{smps_c}')
        clr.append('mediumpurple')
    box_plot(y,leg, clr, append=append)

def cor_box_call(data, slct, append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive box plot based on 
    the selected columns. Comparing different corr mode
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
    columns = data.columns.to_numpy()
    for ccn_c in list(slct.keys()):
        for col in [c for c in columns if f'{ccn_c}_' in c]:
            cols.append(col)
        cols.append(slct[ccn_c])
    # input(cols)
    data = data.dropna(subset = cols)
    for ccn_c in list(slct.keys()):
        smps_c = slct[ccn_c]
        for col in [c for c in columns if (f'{ccn_c}_' in c) |(c ==ccn_c)]:
            y.append(data[col].to_numpy())
            leg.append(f'{col.replace('(cm-3)_cor_setpt', ': ss%=').replace('_','\n')}')
            clr.append('skyblue')
        y.append(data[smps_c].to_numpy())
        leg.append(f'N{smps_c}')
        clr.append('mediumpurple')
    box_plot(y,leg, clr, append=append)

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

def chem_line_plot(x,y,legs, append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    twx = ax.twinx()
    lines = []
    for i in range(len(y)):
        if i == 0:
            L, = ax.plot(x,y[i], label = legs[i])
            lines.append(L)
        else:
            print(legs[i])
            print(y[i])
            print(twx.plot(x,y[i], label = legs[i]))
            L, = twx.plot(x,y[i], label = legs[i])
            lines.append(L)
    labs = [l.get_label() for l in lines]
    leg = ax.legend(lines, labs, loc=0)
    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        if event.dblclick: #Double click to remove all other lines from plot
            legline = event.artist
            for legL in [l for l in list(lined.keys()) if l != legline]:
                origline = lined[legL]
                vis = not origline.get_visible()
                origline.set_visible(vis)
                # Change the alpha on the line in the legend so we can see what lines
                # have been toggled
                if vis:
                    legL.set_alpha(1.0)
                else:
                    legL.set_alpha(0.2)
        else: 
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
    ax.set_ylabel('Concentration [%]')
    twx.set_ylabel('Deviation [%]')
    ax.set_xlabel('Date')
    if append == 0:
        ax.set_title(f"Deviation between CCN and SMPS counts in comparison to {legs[0]}")
    else:
        ax.set_title(f"Deviation between CCN and SMPS counts in comparison to {legs[0]} for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

def chem_scat_plot(x,y,fit_x,fit_y,legs,chem= 'org',append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    for i in range(len(y)):
        L, = ax.plot(x[i],y[i], label = legs[i],ls = '', marker = '*')
        lines.append(L)
    fit_lines = []
    leg = ax.legend()
    # for i in range(len(fit_y)):
    #     fit, = ax.plot(fit_x[i],fit_y[i])
    #     fit_lines.append(fit)
    lined = dict()
    fitlined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline
    # for legline, fitline in zip(leg.get_lines(), fit_lines):
    #     legline.set_picker(5)  # 5 pts tolerance
    #     fitlined[legline] = fitline

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        if event.dblclick: #Double click to remove all other lines from plot
            legline = event.artist
            for legL in [l for l in list(lined.keys()) if l != legline]:
                origline = lined[legL]
                vis = not origline.get_visible()
                origline.set_visible(vis)
                # Change the alpha on the line in the legend so we can see what lines
                # have been toggled
                if vis:
                    legL.set_alpha(1.0)
                else:
                    legL.set_alpha(0.2)
        legline = event.artist
        origline = lined[legline]
        # fitline = fitlined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # fitline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel('Deviation [%]')
    ax.set_xlabel('Composition [%]')
    if append ==0:
        ax.set_title(f"Comparison of CCN and SMPS deviation to {chem} compostion")
    else:
        ax.set_title(f"Comparison of CCN and SMPS deviation to {chem} compostion for {append}")
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
        ax.set_title(f"Comparison of CCN and SMPS data")
    else: 
        ax.set_title(f"Comparison of CCN and SMPS data for {append}")
    ax.set_ylabel('N [#/cm^3]')
    input('Press enter to exit plot...')
    plt.ioff()