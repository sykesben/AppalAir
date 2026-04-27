"""
Date: 2/15/2026
Author: Ben Sykes
Purpose: Process through raw CCN data and remove large spikes and preform basic QA procedures
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import os
from os.path import expanduser 
import datetime as dt 
import scipy as sp
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from CCN_EBAS_convert import ebas_genfile
from pathlib import Path
from scipy.signal import find_peaks
from scipy import stats
from CCN_process import *

def findSpikes(data,col, thresh = 10, width_max = 10):
    vals = data[col].values
    height = np.nanmean(vals)+np.nanmean(vals)*thresh/100
    peaks, _ = find_peaks(vals, height=height, width=[0, width_max])
    bad_indexs = data.index.to_numpy()[peaks]
    data = data.drop(index = bad_indexs)
    return data, peaks

def findZSpikes(data, col, thresh =2):
    vals = data[col].values
    z = np.abs(stats.zscore(data[col],nan_policy='omit'))
    peaks = np.where(z > thresh)[0]
    bad_indices = data.index.to_numpy()[peaks]
    data = data.drop(index =bad_indices)
    return data, peaks

def findZSpikesRoll(data, col, thresh =2, window =20):
    col_mean = data[col].rolling(window=window).mean()
    col_std = data[col].rolling(window=window).std()
    z = (data[col] - col_mean)/col_std
    # input(z)
    peaks = np.where(z > thresh)[0]
    bad_indices = data.index.to_numpy()[peaks]
    data = data.drop(index =bad_indices)
    return data, peaks

def line_plot(x,y,legs):
    plt.ion()
    bad_indexes = []
    fig, ax = plt.subplots(sharex=False)
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
    
    # def onselect(eclick, erelease):
  
    #     # Obtain (xmin, xmax, ymin, ymax) values
    #     # for rectangle selector box using extent attribute.
    #     x1, y1 = eclick.xdata, eclick.ydata
    #     x2, y2 = erelease.xdata, erelease.ydata
    #     print("Date_range: ", f'{y1}-{y2}')
        
    #     date_start = y1
    #     date_end = y2


    fig.canvas.mpl_connect('pick_event', onpick)
    # rect_selector = RectangleSelector(ax, onselect, useblit=True, button=[1], minspanx=5, minspany=5, spancoords='data')
    ax.set_ylabel('N [#/cm^3]')
    ax.set_xlabel('Date')
    ax.set_title(f"Cleaned CCN data")
    input('Press enter to exit plot...')
    plt.ioff()

file_in = r"C:\Users\bensy\Documents\Research\CCN\app_20260101.csv"
# start_file = r"C:\Users\bensy\Documents\Research\CCN_Clean_2026_1min.csv"
file_out = r"C:\Users\bensy\Documents\Research\CCN_Clean_2026_1min.csv"

end_data, rename = readin(file_in)
# print('end')
# start_data,rename = readin(start_file)
# print('start')

data,rename = readin(file_in)

data = data.drop(index= '1970-01-01 00:00:00')
data_clean,peaks = findZSpikesRoll(data, 'N(cm-3)',3,100)
data['check'] = data['N(cm-3)'].to_numpy()
data['check'].iloc[peaks] = np.nan
data.loc[data['N(cm-3)']>5000] = np.nan

data.loc[pd.to_datetime('2026-01-09 17:00:00'): pd.to_datetime('2025-01-09 19:00:00')] = np.nan
data.loc[pd.to_datetime('2026-01-20 17:00:00'): pd.to_datetime('2025-01-20 19:05:00')] = np.nan
data.loc[pd.to_datetime('2026-02-10 13:30:00'): pd.to_datetime('2025-02-10 14:00:00')] = np.nan
data.loc[pd.to_datetime('2026-03-25 15:48:00'): pd.to_datetime('2025-03-25 17:55:00')] = np.nan

date = data.index.to_numpy()
y = [data['N(cm-3)'].to_numpy(),data['check'].to_numpy()]
leg = ['Raw Data', 'Cleaned Data']
line_plot(date,y,leg)

data.iloc[peaks] = np.nan

out = data
input(out)
out.to_csv(file_out)
# if __name__ == '__main__':
