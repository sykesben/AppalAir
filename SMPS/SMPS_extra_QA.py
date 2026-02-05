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


def FindOutliersCOV(data, name, avg_mult = 0.4,size = 10):
    """
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data and the name of a collumn in that dataframe,
    identifies outliers within the named column of the data set via coefficient of variation measurements,
    and returns a dataframe only containing the outliers

    ----------
    Paramaters
    ++++++++++
    data : [Pandas DataFrame] SMPS dataframe with metadata removed
    name : [str] name of the column used to id outliers
    avg_mult : [float] value for deviation check
    size : [float] size of window for rolling operation

    Returns
    ++++++++++
    Outliers : [Pandas DataFrame] outliers compiled into a data frame
    """

    def VarCheck(x):
        if not(x.any ==0):
            avg = np.mean(x)
            if abs(x[-1] -avg) > avg_mult*avg:
                return x[-1]
            else:
                return np.nan
        else:
            return np.nan
    def KeptCheck(x):
        return x.isnull().sum() * 100 / len(x)

    #itterating through each row in the provided data frame
    Outliers = data[name].rolling(window=size).apply(VarCheck,raw =True)
    perc_kept = Outliers.rolling(window='1h').apply(KeptCheck,raw =False)
    data['outliers'] = Outliers
    data['%kept']= perc_kept
    return data   #return the outliers as a series the same size as data[name] with only outliers 

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



filepath = Path(input("\nEnter full path of file you would like to quality assure.\n")) 
data = pd.read_csv(filepath)
print(data.columns)

col1 = 'Aerosol Humidity (%)' #Column names the program is looking at to quality assure
#col2 = 'Geo. Mean (nm)'

data = RemoveOutliers(data,col1)
#data = RemoveOutliers(data,col2)

ParentPath = filepath.parent
name = input('\nEnter the desired name of your combined file and include the file type .csv:\n' \
                     '(This will place the quality assured file in the same folder that held the original file)\n')
data.to_csv(ParentPath / name)