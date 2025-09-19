#%%
#Aydan Gibbs
#Updated 12/31/24
#plots multiple histograms of SMPS Data on the same graph with a legend that displays the days

#Next steps
#-pull files from a folder instead of individually
#-incorporate a for loop to create the graph instead of graphing individually
#-look at data from the before and after the wind turbine went up


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import smps

def main():

    #pulls in file of name f and reads in number of meta lines
    f = r'C:\Users\aydan\VS Code\AppalAIR Code\Raw Data\SMPS_3082002329005_20240611.csv' #test file 1
    f1 = r'C:\Users\aydan\VS Code\AppalAIR Code\Raw Data\SMPS_3082002329005_20240612.csv' #test file 2
    metaDataLines = _get_linecount(f, keyword = 'Scan Number')  #returns the linecount of metadata

    #saves meta data into list meta
    meta = pd.read_table(   #reads in the metadata into a df 'meta'
        f, 
        nrows=metaDataLines, 
        delimiter=',', 
        header=None, 
        encoding='ISO-8859-1',
        on_bad_lines='warn',
        index_col = 0
    ).T.iloc[0,:].to_dict()

    #works with test file 1 to grab the necessary data to display, dataConc and bins
    dataRaw = pd.read_table(                                                                        #reads in the data, skipping over the metaDataLines
        f,
        skiprows = metaDataLines,
        delimiter = ','
    )
    dataConc = dataRaw.get(dataRaw.columns[range(42, 156)])                                         #gets the concentration data from the raw data. assumes 114 columns of concentration data
    SMPSsizeBins = list(dataConc.columns)                                                           #returns a list of strings of the sizebins, assumes 114 columns of data,
    SMPSsizeBins = list(map(float,SMPSsizeBins))                                                    #returns sizebins as a list of floats
    bins = smps.utils.make_bins(lb = 13.3, ub = 805.8,midpoints = np.array(SMPSsizeBins))           #returns a 3xn matrix of size bins with lower, upper and midpoints
    
    #works with test file 2 to, grab the necessary data to display, dataConc1 and bins
    dataRaw1 = pd.read_table(                                                                       #reads in the data, skipping over the metaDataLines
        f1,
        skiprows=metaDataLines,
        delimiter = ','
    )
    dataConc1 = dataRaw1.get(dataRaw1.columns[range(42, 156)])                                      #gets the concentration data from the raw data. assumes 114 columns of concentration data

    #Works with the display data to format both dataConc and dataConc1 on a single graph with a legend
    dates = ["6-11-24", "6-12-24"]                                                                      #creates the dates for the legend
    ax = None                                                                                           #creates the axis to be plotted on
    plot_kws = dict(alpha=0.65, color='r', linewidth=0.)                                                #plot keywords for dataConc 
    ax = smps.plots.histplot(dataConc, bins, ax=ax, plot_kws=plot_kws, fig_kws=dict(figsize=(12, 6)))   #plots dataConc
    plot_kws = dict(alpha=0.65, color='b', linewidth=0.)                                                #plot keywords for dataConc1
    ax = smps.plots.histplot(dataConc1, bins, ax=ax, plot_kws=plot_kws, fig_kws=dict(figsize=(12, 6)))  #plots dataConc1
    ax.legend(dates, loc='best')                                                                        #creates the legend for the graph
    ax.set_ylabel("$dN/dlogD_p \; [cm^{-3}]$")                                                          #lables the y axis
    plt.title('data from 6-11 and 6-12')                                                                #title for graph
    ax.grid()                                                                                           #shows the gridlines
    plt.xlabel('Dp (nm)')                                                                               #lables the x axis
    plt.show()                                                                                          #displays the graph


#returns the line count of the meta data
def _get_linecount(fpath, keyword, delimiter=',', encoding='ISO-8859-1'):
    """
    Return the line number in a file where the first item is 
    `keyword`. If there is no such line, it will return the total 
    number of lines in the file.
    
    :param fpath: A path or url for the file (if a url it must 
        include `http`, and if a file path it must not 
        contain `http`).
    :type fpath: string
    :param keyword: The string to look for.
    :type keyword: string
    :param delimiter: The delimiter between items in the file, 
        defaults to ','.
    :type delimiter: string
    :param encoding: The encoding for the file, defaults to 
        'ISO-8859-1'.
    :type encoding: string
    :return: The line number in the file where the first item is 
        `keyword`.
    :rtype: int
    """
    linecount = 0

    if 'http' in fpath:
        req = requests.get(fpath)

        for line in req.iter_lines():
            startswith = line.decode(encoding).split(delimiter)[0]

            if startswith == keyword:
                break

            linecount += 1
    else:
        with open(fpath, 'r', encoding=encoding) as f:
            for line in f:
                startswith = line.split(delimiter)[0]

                if startswith == keyword:
                    break

                linecount += 1

    return linecount

if __name__:
    main()