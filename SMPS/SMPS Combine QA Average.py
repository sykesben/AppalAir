#Aydan Gibbs
#8/28/25
#Compiles SMPS AIM files in a designated folder and Outputs them into one csv file with no AIM meta data
#Quality assures data and outputs one csv file with no AIM meta data
#Averages data over various time steps and outputs one csv file with no AIM meta data

#Next Steps
#export with metadata
#allow combining to a file that already is without metadata (using a previously combined file)
#include try catch for if there is no outliers in data



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def main():

    #ask user if they would like to combine files
    CombineFilesYN = input("\nWould you like to combine files? (Y/N)\n")
    if CombineFilesYN == 'Y':                                                   #if yes, run CombineFiles
        CombinedDF = CombineFiles()

    #ask user if they would like to QA a file
    QualityAssureFileYN = input('\nWould you like to quality assure a file? (Y/N)\n')
    if QualityAssureFileYN == 'Y':                                              #if yes:
        if CombineFilesYN == 'Y':                                               #and if they combined files previously
                                                                                #ask if user would like to use the combined file
            QualityAssureCombinedFileYN = input('\nWould you like to use the file you just made? (Y/N)\n')
            if QualityAssureCombinedFileYN == 'Y':                              #if yes:
                QADF = QualityAssureFile(CombinedDF, 'N')                       #run QualityAssureFile with previously combined data
            else:                                                               #if no:
                                                                                #prompt for the full path of data they would like to quality assure
                filepath = Path(input("\nEnter full path of file you would like to quality assure.\n"))
                QADF = QualityAssureFile(pd.DataFrame(), filepath)              #run QualityAssureFile with an empty dataframe and the path of the file specified
        else:                                                                   #if they didnt combine a file previously:
                                                                                #prompt for the full path of data they would like to quality assure
            filepath = Path(input("\nEnter full path of file you would like to quality assure.\n"))
            QADF = QualityAssureFile(pd.DataFrame(), filepath)                  #run QualityAssureFile with an empty dataframe and the path of the file specified

    #ask user if they would like to average a file
    AverageFileYN = input("\nWould you like to average a file? (Y/N)\n")
    if AverageFileYN == 'Y':                                                    #if yes:
        if CombineFilesYN == 'Y' or QualityAssureFileYN == 'Y':                 #and if they combined files or QA a file previously
                                                                                #ask if they would like to use the previous file
            AverageCombinedorQualityAssuredFileYN = input("\nWould you like to use the file you just made? (Y/N)\n")
            if AverageCombinedorQualityAssuredFileYN == 'Y' and QualityAssureFileYN == 'Y':
                AverageFile(QADF,'N')                                           #if yes and they used QualityAssureFile:
                                                                                #run AverageFile with QADF
            elif AverageCombinedorQualityAssuredFileYN == 'Y' and CombineFilesYN == 'Y':
                AverageFile(CombinedDF,'N')                                     #if yes and they only combined a file:
                                                                                #run AverageFile with CombinedDF

            else:                                                               #if they dont want to use the file they previously used:
                                                                                #promt for the full path of data they would like to average
                filepath = Path(input("\nEnter full path of file you would like to average.\n"))
                AverageFile(pd.DataFrame(),filepath)                            #run AverageFile with an empty dataframe and the path of file specified
        else:                                                                   #if they didnt use either previous function:
                                                                                #prompt for the full path of data they would like to avearge
                filepath = Path(input("\nEnter full path of file you would like to average.\n"))
                AverageFile(pd.DataFrame(),filepath)                            #run AverageFile with an empty dataframe and the path of file specified




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


    with open(fpath, 'r', encoding=encoding) as f:
        for line in f:
            startswith = line.split(delimiter)[0]

            if startswith == keyword:
                break

            linecount += 1

    return linecount

def CombineFiles():
    """
    8/20/25
    Outputs a dataframe of combined TSI EC 3082 and CPC 3750 SMPS data files within a user specified folder

    :return: A data frame of combined data from the user specified folder
    :rtype: dataframe 
    """

    #create a dataframe to store combined data and metadata
    dataTotal = pd.DataFrame()
    metaTotal = pd.DataFrame() 

    #Get the path to the data folder
    folderpath = Path(input("\nInput the full path of the folder youd like to access:\n"))
    ParentPath = folderpath.parent

    #itterates through each file in the user specified folder and appends them to dataTotal
    for entry in folderpath.iterdir():                                          #looks at each item in the folder
        print(entry)                                                            #prints each file name
        metaDataLines = _get_linecount(entry, keyword = 'Scan Number')          #returns the linecount of metadata

        meta = pd.read_table(                                                   #reads in the metadata into a df 'meta'
                        entry, 
                        nrows=metaDataLines, 
                        delimiter=',', 
                        header=None, 
                        encoding='ISO-8859-1',
                        on_bad_lines='warn',
                        index_col = 0
                    ).T.iloc[0,:].to_dict()

        dataRaw = pd.read_table(                                                #reads in the data, skipping over the metaDataLines
                        entry,
                        skiprows=metaDataLines,
                        delimiter = ','
                    )
        
        dataTotal = dataTotal._append(dataRaw, ignore_index = True)             #append each file to dataTotal
        
        metaTotal = metaTotal._append(meta, ignore_index = True)                #appends each metadata to metaTotal 9/19

    #Convert the "DateTime Sample Start" column to a datetime object
    dataTotal["DateTime Sample Start"] = pd.to_datetime(dataTotal["DateTime Sample Start"], format = 'mixed', dayfirst=True)
    dataTotal = dataTotal.set_index('DateTime Sample Start')                    #now use the datetime object as the new index, this sorts the data by date
    dataTotal = dataTotal.sort_index()

    #saves user data apon request
    CreateFileYN = input('\nWould you like to save this combined file? (Y/N)\n')#promt user to save combined file
    if CreateFileYN == 'Y':                                                     #if yes, create the file at user speified location 
                                                                                #just outside the folder the user specified earlier
        name = input('\nEnter the desired name of your combined file and include the file type .csv:\n' \
                     '(This will place the file just outside the folder you indicated previously with the name you specify)\n')
        dataTotal.to_csv(ParentPath / name)
    print(dataTotal)                                                            #print the combined file for a check

    CreateMetaYN = input('\nWould you like to save this combined metadata? (Y/N)\n')
    if CreateMetaYN == 'Y':                                                     #if yes, create the file at user speified location 
                                                                                #just outside the folder the user specified earlier
        name = input('\nEnter the desired name of your metadata file and include the file type .csv:\n' \
                     '(This will place the file just outside the folder you indicated previously with the name you specify)\n')
        metaTotal.to_csv(ParentPath / name)
    print(metaTotal)                                                            #print the metadata file for a check
    return dataTotal                                                           #return the combined file

def AverageFile(DataDF,FilePath):
    """
    8/28/25
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data with metadata removed.
    Outputs a dataframe of TSI EC 3082 and CPC 3750 SMPS data averaged over a user specified time step

    :param DataDF: TSI EC 3082 and CPC 3750 SMPS dataframe with metadata removed
    :type DataDF: dataframe

    :param filepath: Optional, path of a .csv file that contains TSI EC 3082 and CPC 3750 SMPS
    data with metadata removed, defaults to 'N'
    :type filepath: string

    :return: DataDF averaged over a user specified time step
    :rtype: dataframe 
    """


    if FilePath == "N":                                                         #if there is no filepath use DataDF as dataRaw
        dataRaw = DataDF

        #Takes the statistics, raw, and corrected data columns from the data to then be averaged
        #we do this so that you arnt trying to average N/A data, or text data
        AllColumns = list(dataRaw.columns)                                      #lists all the columns
        StatsHeaders = AllColumns[33:40] + AllColumns[41:500]                   #selects statistics columns, raw, and corrected data from the list
        print(StatsHeaders)
        dataRaw =  dataRaw[StatsHeaders]                                        #uses just the previously selected columns
        dataRaw.index = pd.to_datetime(dataRaw.index)                           #turn the index back into a date time object (this was undone some how previously)

    else:                                                                       #if there is a filepath, load that file as dataRaw
        if FilePath.is_file():
                        print(FilePath)
                        dataRaw = pd.read_csv(FilePath)

                        #Convert the "DateTime Sample Start" column to a datetime object
                        dataRaw["DateTime Sample Start"] = pd.to_datetime(dataRaw["DateTime Sample Start"], format = 'mixed', dayfirst=True)

                        #Takes the statistics, raw, and corrected data columns from the data to then be averaged
                        #we do this so that you arnt trying to average N/A data, or text data
                        AllColumns = list(dataRaw.columns)                      #lists all the columns
                        StatsHeaders = AllColumns[33:40] + AllColumns[41:500]   #selects statistics columns, raw, and corrected data from the list
                        StatsHeaders.append('DateTime Sample Start')            #add the time stamps to the list
                        dataRaw =  dataRaw[StatsHeaders]                        #uses just the previously selected columns   
                        dataRaw = dataRaw.set_index('DateTime Sample Start')    #now use the datetime object as the new index, this sorts the data by date

    #averages the data based on a user inputted time step                     
    StepSize = input('\nEnter a time step for the averaging.\n'                 #user inputs the desired time step
                     'To format the time step include a number followed by the unit of time, ex. 5h = 5 hours time step\n'
                     's = seconds, min = minutes, h = hours, d = days, W = weeks, M = months\n')
    dataRaw = dataRaw.resample(StepSize).mean()                                 #averages the data over the designated time step
    print(dataRaw)                                                              #displays data so you can check its the timestep you wanted

    #saves user data apon request
    CreateFileYN = input('\nWould you like to save this averaged file? (Y/N)\n')#prompt the user to save the averaged file
    if CreateFileYN == 'Y':                                                     #if yes: create a file
        if FilePath == 'N':                                                     #if there was not a filepath passed at the beginning, 
                                                                                #prompt for a filepath to save the data to
            name = input('\nEnter the full path for your averaged file and include the file type .csv:\n')    
            dataRaw.to_csv(name)
        else:                                                                   #if there was a filepath passed
                                                                                #prompt for a name for the file, and save in the same folder as the filepath
            ParentPath = FilePath.parent
            name = input('\nEnter the desired name for your averaged file and include the file type .csv:\n'
                         '(This will place the file in the same folder as the file you just averaged with the name you specify)\n')    
            dataRaw.to_csv(ParentPath / name)


def QualityAssureFile(DataDF ,filepath = 'N'):
    """
    8/20/25
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data with metadata removed.
    Outputs a dataframe of TSI EC 3082 and CPC 3750 SMPS data with outliers removed 
    by FindOutliersAverage and FindOutliersRange removed.

    :param DataDF: TSI EC 3082 and CPC 3750 SMPS dataframe with metadata removed
    :type DataDF: dataframe

    :param filepath: Optional, path of a .csv file that contains TSI EC 3082 and CPC 3750 SMPS
    data with metadata removed, defaults to 'N'
    :type filepath: string

    :return: DataDF with outliers removed
    :rtype: dataframe 
    """
    
    
    if filepath == 'N':                                                         #if there is no filepath use DataDF as dataRaw
        dataRaw = DataDF
    else:                                                                       #if there is a filepath, load that file as dataRaw
        if filepath.is_file():
            print('\nFile to be quality assured: ' + str(filepath))
            dataRaw = pd.read_table(
                filepath,
                delimiter = ',')

            #Convert the "DateTime Sample Start" column to a datetime object
            dataRaw["DateTime Sample Start"] = pd.to_datetime(dataRaw["DateTime Sample Start"], format = 'mixed', dayfirst=True)
            dataRaw = dataRaw.set_index('DateTime Sample Start')                #now use the datetime object as the new index, this sorts the data by date


    #use outlier functions to remove outliers
    HumidityRangeOutliers = FindOutliersRange(dataRaw, 'Aerosol Humidity (%)', 0, 40)
    dataRaw = RemoveOutliers(dataRaw, HumidityRangeOutliers)
    AverageOutliers = FindOutliersAverage(dataRaw, 'Geo. Mean (nm)')
    dataRaw = RemoveOutliers(dataRaw, AverageOutliers)

    #use simple pandas functions to remove data that are not 'Normal Scans'
    dataRaw = dataRaw[dataRaw['Detector Status'] == 'Normal Scan']
    dataRaw = dataRaw[dataRaw['Classifier Errors'] == 'Normal Scan']

    #saves data upon request
    CreateFileYN = input('\nWould you like to save this QA file? (Y/N)\n')      #prompt the user to save the QA file
    if CreateFileYN == 'Y':                                                     #if yes: create a file
        if filepath == 'N':                                                     #if there was not a filepath passed at the beginning, 
                                                                                #prompt for a filepath to save the data to
            name = input('\nEnter the full path for your QA file and include the file type .csv:\n')    
            dataRaw.to_csv(name)
        else:                                                                   #if there was a filepath passed
                                                                                #prompt for a name for the file, and save in the same folder as the filepath
            ParentPath = filepath.parent
            name = input('\nEnter the desired name of your QA file and include the file type .csv:\n' \
                         '(This will place the file same folder as the file you just QA with the name you specify)\n')
            dataRaw.to_csv(ParentPath / name)
    print(dataRaw)
    return dataRaw                                                              #return QA dataframe


def FindOutliersAverage(Data, DataType):
    """
    8/19/25
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data and the name of a collumn in that dataframe,
    identifies outliers within the named collumn of the data set via averaging surrounding data and comparing,
    and returns a dataframe containing the outliers

    :param Data: TSI EC 3082 and CPC 3750 SMPS dataframe with metadata removed
    :type Data: dataframe

    :param dataType: name of the column used to id outliers
    :type dataType: string

    :return: outliers compiled into a data frame
    :rtype: dataframe 
    """
    
    #define local variables

    #rows to be avearaged together
    row1 = 0
    row2 = 0
    row3 = 0
    row4 = 0
    row5 = 0
    row6 = 0
    row7 = 0
    row8 = 0
    row9 = 0

    AveMultiplyer = .4                                                          #the initial % of the average of compared values that is used to select if a value is an outlier or not (can be changed)
    Outliers = pd.DataFrame()                                                   #a blank data frame to store the outlier terms

    #itterating through each row in the provided data frame
    for index, row in Data.iterrows():
        row1 = row[DataType]                                                    #assign the current row to row1 term


        #if all row terms have been populated with a value from the provided data look for an outlier, 
        #this makes sure that when averaging you are averaging over data instead of the preset 0's
        if row1 != 0 and row2 != 0 and row3 != 0.0 and row4 != 0.0 and row5 != 0.0 and row6 != 0 and row7 != 0 and row8 != 0 and row9 != 0:

            Ave = (row2+row3+row4+row5+row6+row7+row8+row9)/8                   #average all the row terms together except row1 (the term checked for being an outlier)

            if abs(row1 - Ave) > AveMultiplyer*Ave:                             #if the difference of the average and the check term (row1) is larger than the product 
                                                                                #of the multiplyer and the average then:
                Outliers = pd.concat([Outliers, pd.DataFrame([row])])           #add the check term (row1) to the outliers dataframe

                                                                                #if you are removing the check term (row1) then it will not be incorperated to
                                                                                #the averaging and only row1 will be updated with the next term in the provided data

            #if the check term (row1) is not an outlier then update all the terms
            #row1 becomes the next term in the provided data, row2 becomes row1, row3 becomes row2 and so on
            else:
                row9 = row8
                row8 = row7
                row7 = row6
                row6 = row5
                row5 = row4
                row4 = row3
                row3 = row2
                row2 = row1
        
        #if all of the rows are not populated with data from the provided data then update all the terms
        #row1 becomes the next term in the provided data, row2 becomes row1, row3 becomes row2 and so on
        else:
            row9 = row8
            row8 = row7
            row7 = row6
            row6 = row5
            row5 = row4
            row4 = row3
            row3 = row2
            row2 = row1

    return Outliers                                                             #return the dataframe of Outliers

def FindOutliersRange(Data,Column,Min,Max):
    """
    8/19/25
    Takes in a data frame of TSI EC 3082 and CPC 3750 SMPS data and identifies outliers 
    within the data set that are outside of the Min and Max indicated values. 
    Returns a dataframe containing the outliers.

    :param Data: TSI EC 3082 and CPC 3750 SMPS data with metadata removed
    :type Data: dataframe
    :param Column: name of the column used to id outliers
    :type Column: string
    :param Max: data greater than this will be an outlier
    :type Max: intiger
    :param Min: data less than this will be an outlier
    :type Min: intiger

    :return: dataframe of Outliers
    :rtype: dataframe
    """

    #creates a dataframe of all data in the designated Column that is greater than Max or less than Min
    Outliers = pd.concat([Data[Data[Column]<Min], Data[Data[Column]>Max]])

    return Outliers                                                             #return the dataframe of outliers

def RemoveOutliers(Data, Outliers):
    """
    8/20/25
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data and an outlier dataframe

    :param Data: TSI EC 3082 and CPC 3750 SMPS data with metadata removed
    :type Data: dataframe
    :param Outliers: Outlier data, returned from either FindOutliersRange or FindOutliersAverage, 
    you would like to remove from Data
    :type Outliers: dataframe

    :return: Original Data with all data in Outliers removed from it
    :rtype: dataframe
    """
    columns = Data.columns                                                      #Save column names for correcting later

    #merge the Outliers with Data and mark where values overlap with flag "left_only"
    #to keep the date time index you must use left_index = True, and right_index = True otherwise it will apply its own intiger index
    #but this will double the columns which are removed later
    Data = Data.merge(Outliers, how='left', indicator = True, left_index = True, right_index = True)
    
    Data = Data[Data['_merge'] == 'left_only']                                  #flag overlapped data with "_merge"
    DataNoOutliers = Data.drop(columns = '_merge')                              #remove all data with the overlap flag "_merge"

    #remove artifacts generated by the .merge function
    DataNoOutliers = DataNoOutliers.drop(columns = DataNoOutliers.columns[(len(columns)):2*len(columns)])
    DataNoOutliers.columns = columns

    return DataNoOutliers                                                       #return Data with outliers removed


#if the program exists, run it
if __name__:
    main()