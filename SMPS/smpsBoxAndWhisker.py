#%% 
#Aydan Gibbs
#1/22/25
#box and wisker plots for smps


#Next Steps
#sorting files by thier names to put them in order on the graph

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import smps
import os
import requests
from pathlib import Path

def main():

    dataTotal = pd.DataFrame()
    dateTotal = pd.DataFrame()
    dataType = input('Enter the type of data you would like to display (Geo. Mean (nm), Geo. Std. Dev, Total Concentration (#/cmÂ³)): ')
    fileNames = []
    csvpath = Path(input("\nInput the full path of the csv youd like to access:\n"))
    
    
    
    #directory = (r'C:\Users\aydan\VS Code\AppalAIR Code\Raw Data\Compiled')                                                 #opens the compiled folder
    for entry in os.listdir(csvpath):                                                                                     #looks at each item in the folder
        file = os.path.join(csvpath, entry)                                                                               #selects one item
        if os.path.isfile(file):                                                                                            #if the item is a file it:
                    print(file)                                                                                             #prints the file name
                    dataJoin = pd.DataFrame()                                                                               #initializes an empty data frame for joining
                    dateJoin = pd.DataFrame()
                    dataRaw = pd.read_table(                                                                                #reads in the data, skipping over the metaDataLines
                        file,
                        delimiter = ','
                    )

                    file = file.split('\\')                                                                                 #formats the name of the file to be just one number
                    file = file[-1]
                    file = file.split('.')
                    file = file[0]
                    file = file.split('_')
                    file = file[0]

                    fileNames.append(file)                                                                                  #adds file names to a list
                    dateJoin[file] = dataRaw['DateTime Sample Start']
                    
                    dataJoin[file] = dataRaw[11:263].sum()                  #was dataJoin[file] = dataRaw[dataType]
                                                                      #puts current file data of dataType into dataJoin data frame
                    
                    dataTotal = dataTotal.join(dataJoin, how = 'outer')                                                     
                    #joins dataTotal with dataJoin, fixes truncation of data to first item issue
                    dateTotal = dateTotal.join(dateJoin, how = 'outer')

    fileNames = fileNames[3:] + fileNames[:3]
    
    dataTotal.plot(x = dataTotal,  y = dateTotal, color='blue', linestyle='--', marker='o', label='My Line')                   
    #changing was dataTotal.boxplot(column = fileNames)
    
    plt.xlabel('Month_Year')
    plt.ylabel(dataType)
    plt.title('App_SMPS' + dataType)
    plt.show()

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