"""
Date: 1/29/2026
Author: Ben Sykes
Purpose: Preforming a rolling check of the Coefiecient of Variation as a way to optimize
SMPS processing.
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from pathlib import Path


def FindOutliersRolling(data, name, name_out='', avg_mult = 0.4,size = 10):
    """
    Takes in a dataframe of TSI 3938 SMPS data and the name of a column in that dataframe to process,
    identifies outliers within the named column of the data set via coefficient of variation measurements,
    and returns the original dataframe with outliers marked

    ----------
    Parameters
    ++++++++++
    data : [Pandas DataFrame] SMPS dataframe with metadata removed
    name : [str] name of the column used to id outliers
    name_out : [str] name of the flag column for outputting (default = '')
    avg_mult : [float] value for deviation check (default = 0.4)
    size : [float] size of window for rolling operation (default = 10)

    Returns
    ++++++++++
    Outliers : [Pandas DataFrame] outliers compiled into a data frame
    """
    # generate a forward mean
    forward_mean = (data[name].shift(-1)        # Shift all rows back by 1
                    .rolling(window=size-1)     # Look forward at the next N-1 rows
                    .mean()                     # Mean these next N-1 rows (mean automatically placed at the index of the last row)
                    .shift(-(size-2)))          # Shift mean back from last index to first index
    # Check if the value of each scan exceeds the mean value of the next N-1 scans by more than a certain amount
    outliers = ((data[name] - forward_mean).abs() >avg_mult * forward_mean) 
    if name_out == '': name_out = name + ' outliers'
    data[name_out] = outliers
    return data   #return dataframe with the outlier added in

def RemoveOutliers(data, name, avg_mult = 0.4,size= 10):
    """
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data and the name of a collumn in that dataframe,
    identifies outliers within the named column of the data set via coefficient of variation measurements,
    and returns a dataframe without the outlier rows

    ----------
    Paramaters
    ++++++++++
    data : [Pandas DataFrame] SMPS dataframe with metadata removed
    name : [str] name of the column used to id outliers
    avg_mult : [float] value for deviation check
    size : [float] size of window for rolling operation

    Returns
    ++++++++++
    data : [Pandas DataFrame] Updated data with outliers removed
    """
    start_check = CheckWindow(data,name)
    def VarCheck(x):
        if not(x.any ==0):
            avg = np.mean(x)
            if abs(x[-1] -avg) > avg_mult*avg:
                return 1
            else:
                return 0
        else:
            return 0
    def KeptCheck(x):
        return np.count_nonzero(x==0)/len(x)*100
    #itterating through each row in the provided data frame
    print(data[name].rolling(window=size).apply(VarCheck,raw =True))
    data['out'] = data[name].rolling(window=size).apply(VarCheck,raw =True)
    data['% kept'] = data['out'].rolling(window='1h').apply(KeptCheck,raw =True)
    data = data[data[name].rolling(window=size).apply(VarCheck,raw =True)==0]
    data = data.drop(columns ='out')
    return data   #return the dataframe with the outlier rows removed
    
def CheckWindow(data, name,start= 0,avg_mult = 0.4,size =10):
    """
    Takes in the first N[size] points in a column to verify the validity of the
    QA check

    ----------
    Paramaters
    ++++++++++
    data : [Pandas DataFrame] SMPS dataframe with metadata removed
    name : [str] name of the column used to id outliers
    avg_mult : [float] value for deviation check (default = 0.4)
    size : [float] size of window for rolling operation (defualt )
    

    Returns
    ++++++++++
    Outliers : [Pandas DataFrame] outliers compiled into a data frame
    """
    total_avg = data[name].mean()
    start_avg = data[name].iloc[start:start+size].mean()
    valid = np.abs(start_avg-total_avg) > avg_mult*total_avg
    return valid