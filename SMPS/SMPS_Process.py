# Ben Sykes built on work by Ethan Parkhurst and Aydan Gibbs
# 6/25/26
# Compiles SMPS AIM files in a designated folder and Outputs them into 2 csv files, one at ambient and one at standard
# Averages data over various time steps and outputs one csv file with no AIM meta data

import pandas as pd
import numpy as np
from pathlib import Path, PurePath
from os.path import expanduser 
from SMPS_EBAS_convert import ebas_genfile

'''If ran as primary user, dev should be set to true, and paths should be specified'''
dev = True
cut_RH = 40
year = '2026'                                                               # year to pull from
if year != '':                                                              # if year is provided, set start and end date to only contain dates within the year
    start = pd.to_datetime(f'01/01/{year}')
    end = pd.to_datetime(f'12/31/{year}')
conc = 'num'                                                                # if set up similar to my set up, either num for number or vol for volume concentration
conc_dict = {'vol':'Volume', 'num':'Number', 'mass': 'Mass'}
inFolderpath = expanduser(f'~\Documents\Research\SMPS_{year}\Raw_{conc}')   # folder that contains all the input files and where the output files will be placed
combo_name = fr'{year}_SMPS{conc}'                                          # the name for the outputted raw combined file, .csv extension not included
QA_name = fr'{year}_SMPS{conc}_QA'                                          # the name for the outputted combined QA file, .csv extension not included
TimeStep = 'h'                                                              # timestep to resample data to
out_name = fr'{year}_SMPS_{conc_dict[conc]}SizeDist_1hr'                    # Finalized output file, .csv extension not included
             #2024                                                          # suspected Bad dates to remove from cleaned processing file
bad_dates = [[pd.to_datetime('6/5/2024 20:02:00'),pd.to_datetime('6/6/2024 00:00:00')],[pd.to_datetime('6/5/2024 20:02:00'),pd.to_datetime('6/7/2024 00:00:00')],
             #2025
             [pd.to_datetime('7/10/2025 14:37:00'),pd.to_datetime('7/10/2025 18:20:00')],[pd.to_datetime('7/11/2025 12:48:00'),pd.to_datetime('7/11/2025 17:02:00')],[pd.to_datetime('7/14/2025 13:16:00'),pd.to_datetime('7/15/2025 19:42:00')],
             [pd.to_datetime('7/16/2025 14:30:00'),pd.to_datetime('7/16/2025 20:47:00')],[pd.to_datetime('8/7/2025 13:33:00'),pd.to_datetime('8/7/2025 18:18:00')],[pd.to_datetime('8/7/2025 13:33:00'),pd.to_datetime('8/7/2025 18:18:00')],
             [pd.to_datetime('12/4/2025 17:40:00'),pd.to_datetime('12/6/2025 13:09:00')],
             #2026
             [pd.to_datetime('1/19/2026 10:47:00'),pd.to_datetime('1/19/2026 17:47:00')],[pd.to_datetime('3/20/2026 15:20:00'),pd.to_datetime('3/20/2026 17:37:00')],[pd.to_datetime('04/17/2026 00:00:00'),pd.to_datetime('06/18/2026 16:00:00')],]

def main():
    '''====1 Process Files 1====+++
    Read in the minute averaged files and process them to completeness. If ran as primary user,
    the process is simplified. If dev is false, filepaths can be specified during run time. 
    General Process:
     I : CombineFiles() <- no inputs, pulls all raw files from specified folder and returns combined data
     II : QualityAssureFile(CombinedDF) <- combined data as input, returns quality assured file
     III : AverageFile(QADF) <- Quality assured data as input, returns raw averaged data, cleaned averaged data,
                                stp averaged data, and data with all necessary data for EBAS NASA ames file
     IV : SMPS_EBAS(EBAS, EBAS_Folder) <- takes in EBAS datafile and folder path and outputs converted NASA ames file to path
    +++====1 Process Files 1===='''
    global dev
    if not dev:
       dev_choice = input("\nHave you specified filepaths at the top of the code? (Y/N)\n")
       if dev_choice == 'Y': dev = True
    if dev:                                                                                                 # if ran as a primary user
        CombinedDF = CombineFiles()
        QADF = QualityAssureFile(CombinedDF)
        raw, clean, stp, EBAS = AverageFile(QADF)
        if conc =='num':
            SMPS_EBAS(EBAS, Path(inFolderpath).with_name('ACTRIS'))

    if not dev:
        while True: 
            task = input('What task are you hoping to preform?\n[opts: Combine, QA, Average, STP, All]')
            if task.lower() == 'all':                                                                       # perform all functions on SMPS data
                CombinedDF = CombineFiles()
                QADF = QualityAssureFile(CombinedDF)
                AverageFile(QADF,'N')
            elif task.lower() == 'combine':                                                                 # ask user if they would like to combine files if yes, run CombineFiles
                CombinedDF = CombineFiles()
            elif task.lower() == 'qa':
                QAfilepath = Path(input("\nEnter full path of file you would like to quality assure.\n"))   # prompt for the full path of data they would like to quality assure
                QADF = QualityAssureFile(QAfilepath)                                                        # run QualityAssureFile with an empty dataframe and the path of the file specified
            elif task.lower() == 'average':
                Avgfilepath = Path(input("\nEnter full path of file you would like to average.\n"))         # prompt for the full path of data they would like to average
                AverageFile(pd.DataFrame(),Avgfilepath)                                                     #   run AverageFile with an empty dataframe and the path of file specified
            elif task.lower() == 'stp':
                STPfilepath = Path(input("\nEnter full path of file you would like to convert to STP.\n"))
                Data = ConvertToSTP(STPfilepath)
            else: 
                print(f'{task} not found in list [Combine, QA, Average, STP, All]')
                continue #jump back to top of list

def _get_linecount(fpath, keyword, delimiter=',', encoding='ISO-8859-1'):
    """
    Return the line number in a file where the first item is 
    `keyword`. If there is no such line, it will return the total 
    number of lines in the file.
    ----------
    Parameters
    ++++++++++
    fpath : [str] A path or url for the file (urls must include `http`, paths must not)
    keyword : [str] The string to look for.
    delimiter : [str] The delimiter between items in the file(default = ',')
    encoding : [str] The encoding for the file (default = 'ISO-8859-1')

    Returns
    +++++++
    return : [int] The line number in the file where the first item = `keyword`
    """
    linecount = 0
    with open(fpath, 'r', encoding=encoding) as f:
        for line in f:
            startswith = line.split(delimiter)[0]

            if startswith == keyword:
                break

            linecount += 1
            if linecount > 53:
                linecount = 0
                break
    print(linecount)
    return linecount

def CombineFiles():
    """
    Outputs a dataframe of combined TSI EC 3082 and CPC 3750 SMPS data files within a user specified folder
    ----------
    Returns
    ++++++++++
    return : [Dataframe] A data frame of combined data from the user specified folder
    """
    #create a dataframe to store combined data
    dataTotal = pd.DataFrame()
    metaTotal = pd.DataFrame()

    #Get the path to the data folder
    if dev:
        folderpath = Path(inFolderpath)
    else:
        folderpath = Path(input("\nInput the full path of the folder you'd like to access:\n"))
    ParentPath = folderpath.parent

    #itterates through each file in the user specified folder and appends them to dataTotal
    for entry in folderpath.iterdir():                                          #looks at each item in the folder
        print(entry)                                                            #prints each file name
        if "METADATA" in str(entry):
            meta = pd.read_table(                                                   #reads in the metadata into a df 'meta'
                    entry, 
                    #nrows=metaDataLines, 
                    delimiter=',', 
                    #header=None, 
                    #encoding='ISO-8859-1',
                    #on_bad_lines='warn',
                    index_col = 0)
        else:
            metaDataLines = _get_linecount(entry, keyword = 'Scan Number' or 'DateTime Sample Start')          #returns the linecount of metadata
            meta = pd.read_table(entry,nrows=metaDataLines,delimiter=',',header=None,
                                 encoding='ISO-8859-1',on_bad_lines='warn',index_col = 0
                                ).T.iloc[0,:]  #reads in the metadata into a df 'meta'
            dataRaw = pd.read_table(entry,skiprows=metaDataLines,delimiter = ',')#reads in the data, skipping over the metaDataLines
            dataRaw['Units'] = meta['Units']
            dataRaw['MCC'] = 'MCC'*int(bool(meta['Multiple Charge Correction'].capitalize()))
            dataRaw['DLC'] = 'DLC'*int(bool(meta['Diffusion Loss Correction'].capitalize()))
            dataRaw['Corrections'] = dataRaw['MCC'] +','+ dataRaw['DLC']
            dataRaw.drop(columns=['MCC','DLC'], inplace=True)
        dataTotal = dataTotal._append(dataRaw, ignore_index = True)             #append each file to dataTotal
        metaTotal = metaTotal._append(meta, ignore_index = True)                #appends each metadata to metaTotal 

    #Convert the "DateTime Sample Start" column to a datetime object
    dataTotal['DateTime Sample Start'] = pd.to_datetime(dataTotal["DateTime Sample Start"], format = 'mixed', dayfirst=True)
    dataTotal = dataTotal.set_index('DateTime Sample Start')                    # now use the datetime object as the new index, this sorts the data by date
    dataTotal.index = dataTotal.index.rename('Datetime(UTC)')                   # rename index to be identical to CCN datafiles for easier merging
    dataTotal = dataTotal.sort_index()
    nFiles = len(list(folderpath.iterdir()))                                    # number of files combined
    print(f'{nFiles} files combined for period {np.datetime_as_string(dataTotal.index.to_numpy()[0], unit='h')} - {np.datetime_as_string(dataTotal.index.to_numpy()[-1], unit='h')}')
    if (start != '')&(end!=''):                                                 # ensure data comes from specified year
        mask = (dataTotal.index >= start) & (dataTotal.index <= end)                
        dataTotal = dataTotal.loc[mask]
    #saves user data upon request
    if dev and (combo_name!=''):  
        final_m_name = combo_name + 'METADATA.csv'
        final_d_name = combo_name + '.csv'
        dataTotal.to_csv(ParentPath / final_d_name)
        metaTotal.to_csv(ParentPath / final_m_name) 
        print(f'2 files generated at {ParentPath}')
    elif not dev:
        name = input('\nEnter the desired name of your combined file. [DO NOT INCLUDE .CSV]:\n' \
                    '(This will place the file just outside the input folder)\n')
        final_m_name = name + 'METADATA.csv'
        final_d_name = name + '.csv'
        dataTotal.to_csv(ParentPath / final_d_name)
        metaTotal.to_csv(ParentPath / final_m_name)                         #creates the csv file with name plus META with all the metadata from the combination
        print(f'2 files generated at {ParentPath}')
    return dataTotal                                                            #return the combined file

def QualityAssureFile(Data, folder = ''):
    """
    Takes in a dataframe of TSI EC 3082 and CPC 3750 SMPS data with the metadata header removed.
    Generates flags for the dataframe based on ACTRIS recomendations
    using FindOutliersAverage and FindOutliersRange removed.
    ----------
    Parameters
    ++++++++++
    Data : [DataFrame or path-like] Either dataframe or Path to a .csv file of TSI SMPS with metadata header removed 
    folder : [path-like] Path to folder for outputted data

    Returns
    +++++++
    return : [DataFrame] Data with outliers flagged
    """
    if (isinstance(Data, str))|(isinstance(Data, PurePath)):                                    # If data is a file path, read in Data to dataframe
        print('\nFile to be quality assured: ' + str(Data))
        dataRaw = pd.read_table(Data,delimiter = ',')
        dataRaw['Datetime(UTC)'] = pd.to_datetime(dataRaw['Datetime(UTC)'], format = 'mixed')   # Convert the "Datetime(UTC)" column to a datetime object
        dataRaw = dataRaw.set_index('Datetime(UTC)')                                            # Now use the datetime object as the new index, this sorts the data by date
    else:
        dataRaw = Data                                                                          # Assumed to be a DataFrame, will break otherwise
    # Use outlier functions to flag outliers
    # flag if the humidity in either the Sample or Sheath line is greater than 40 or less than 0
    dataRaw = FindOutliersRange(dataRaw, 'Aerosol Humidity (%)', 0, cut_RH, name_out='Sample RH Flag')
    dataRaw = FindOutliersRange(dataRaw, 'Sheath Relative Humidity (%)', 0, cut_RH, name_out='Sheath RH Flag')
    # flag if the geometric means deviates by more than 40% from the next 9 scans
    dataRaw = FindOutliersRolling(dataRaw, 'Geo. Mean (nm)', name_out='Size Flag')
    if conc =='num': # if number concentration
        # flag if the conc deviates by more than 40% from the next 9 scans or is outside normal values
        dataRaw = FindOutliersRolling(dataRaw, 'Total Concentration (#/cm³)', name_out='N0 Flag')
        dataRaw = FindOutliersRange(dataRaw, 'Total Concentration (#/cm³)', 20, 7000, name_out='N1 Flag')
    elif conc =='vol':
        # flag if the conc deviates by more than 40% from the next 9 scans
        dataRaw = FindOutliersRolling(dataRaw, 'Total Concentration (nm³/cm³)', name_out='V Flag')
    elif conc =='mass':
        # flag if the conc deviates by more than 40% from the next 9 scans
        dataRaw = FindOutliersRolling(dataRaw, 'Total Concentration (µg/m³)', name_out='M Flag')
    # Remove data that are not 'Normal Scans'
    dataRaw = dataRaw[dataRaw['Detector Status'] == 'Normal Scan']
    dataRaw = dataRaw[dataRaw['Classifier Errors'] == 'Normal Scan']

    # Saves data upon request
    if dev and (QA_name!=''): 
        path_out = Path(inFolderpath).with_name(QA_name+'.csv')
        dataRaw.to_csv(path_out)
        print(f'QA file generated at {path_out.parent}')
    elif not dev: 
        if (isinstance(Data, str))|(isinstance(Data, PurePath)):                # If there was a filepath passed, prompt for a name for the file, and save in the same folder as the filepath
            if folder == '': folder = Path(Data).parent
            name = input('\nEnter the desired name of your QA file. [DO NOT INCLUDE .CSV]:\n' \
                        f'(This will place the file within {folder})\n')
            dataRaw.to_csv(Path(Data).with_name(name+'.csv')) 
            print(f'QA file generated at {folder}')
        else: 
            if folder =='':
                path_out = input('\nEnter the full path for your QA file. [DO NOT INCLUDE .CSV]:\n')  
            else: 
                name = input('\nEnter the desired name of your QA file. [DO NOT INCLUDE .CSV]:\n' \
                        f'(This will place the file within {folder})\n')
                path_out = folder+'/'+QA_name+'.csv'
            path_out = Path(path_out)
            dataRaw.to_csv(name)
            print(f'QA file generated at {path_out.parent}')
    return dataRaw                                                                  # Return QA dataframe

def AverageFile(Data,folder =''):
    """
    Takes in a dataframe of TSI SMPS data with metadata removed and averages over a user specified time step
    ----------
    Parameters
    ++++++++++
    Data : [DataFrame or path-like] Either dataframe or Path to a .csv file of TSI SMPS with metadata header removed 
    folder : [path-like] Path to folder for outputted data

    Returns
    +++++++
    return : [DataFrame] Data averaged over a user specified time step
    """
    def quantile(q=0.5, **kwargs):
        def f(series):
            return series.quantile(q, **kwargs)
        return f

    if isinstance(Data, pd.DataFrame):                                          # If data is a DataFrame
        dataRaw = Data
        dataRaw.index = pd.to_datetime(dataRaw.index)                           # ensure the index is a datetime
    if (isinstance(Data, str))|(isinstance(Data, PurePath)):                    # if data is a file path, read in Data to dataframe 
        dataRaw = pd.read_csv(Data)
        dataRaw["Datetime(UTC)"] = pd.to_datetime(dataRaw["Datetime(UTC)"], format = 'mixed')
        dataRaw = dataRaw.set_index('DateTime(UTC)')                            # create and set a datetime index

    cols = dataRaw.columns.to_numpy()
    bad_cols = ['Scan Number','Test Name','Detector Status','Classifier Errors','Communication Status',
                'Neutralizer Status','Detector Inlet Flow (L/min)','Detector Counting Flow (L/min)','Impactor Flow (L/min)',
                'Impactor D50 (nm)','Sheath Flow (L/min)','Scan Direction','HV Polarity','DMA V Ramping Up (TUP) (s)',
                'DMA V Ramping Down (TDOWN) (s)','DMA Column transit time Tf (s)','DMA Exit to Optical Detector Td (s)',
                'DMA at Low Voltage (TLOW) (s)','DMA at High Voltage (THIGH) (s)','Adjustment (s)','Dilution Factor','Aerosol Density (g/cm³)',
                'Reserved 1','Low Voltage (V)','High Voltage (V)','Wide Range Scan Mode',]
    dataRaw = dataRaw.drop(columns=bad_cols,errors='ignore')                    # Drop columns, if columns dont exist, skip 
    dataRaw.columns = dataRaw.columns.str.replace("_", "raw", regex=True)
    # pull out and sort bin size columns 
    numsmps = [s for s in dataRaw.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())&('raw' not in s)]
    numsmps = sorted(numsmps, key=lambda x: float(x))
    meta_cols = [c for c in dataRaw.columns if c not in numsmps]
    reordered_cols = [c for c in dataRaw.columns if c not in numsmps]
    reordered_cols.extend(numsmps) 
    dataRaw = dataRaw[reordered_cols] # reorder columns to put numeric columns after metadata columns
    new_cols = {}
    # rename diameter columns 
    for col in dataRaw.columns.to_numpy():
        if col in numsmps:
            # ensure all bin diameters are 2 decimal places to correctly merge
            float_num = float(col)
            new_cols[col] = f"{float_num:.2f}"
        else:
            new_cols[col] = col
    dataRaw = dataRaw.rename(columns=new_cols)
    # pull out and mean duplicate numeric columns
    df_num = dataRaw[numsmps].copy()
    df_num.columns = df_num.columns.astype(str)
    df_num = df_num.T.groupby(level=0).mean().T
    resort = sorted(df_num.columns.to_numpy(), key=lambda x: float(x))
    df_num = df_num[resort]
    df_meta = dataRaw[meta_cols].copy()
    dataRaw = pd.concat([df_meta, df_num], axis=1)

    #averages the data based on a user inputted time step         
    if (dev & (TimeStep== ''))|(not dev):            
        StepSize = input('\nEnter a time step for the averaging.\n'             # user inputs the desired time step
                        'To format the time step include a number followed by the unit of time, ex. 5h = 5 hours time step\n'
                        's = seconds, min = minutes, h = hours, d = days, W = weeks, M = months\n')
    else: 
        StepSize = TimeStep
    dataClean = FinalCleaning(dataRaw)                                          # clean data by removing flags
    # resample raw data to 'StepSize' timestep
    numeric = dataRaw.select_dtypes(include=np.number).columns
    string = dataRaw.select_dtypes(include=object).columns
    agg_dict = {**{c: 'mean' for c in numeric},**{c: 'first' for c in string}}
    dataRaw = dataRaw.resample(StepSize).agg(agg_dict)                          # averages the numbers and returns the first string over the designated time step
    mask = ((dataRaw == 0) | (dataRaw.isna())).all(axis=1)
    dataRaw = dataRaw.loc[~mask]
    flags = [col for col in dataRaw.columns.to_numpy() if 'Flag' in col]        # round flags to ensure binary representation
    for flag in flags:
        dataRaw[flag] = np.round(dataRaw[flag].to_numpy())
    flags_check = dataRaw[flags].to_numpy()
    flag_code = [int(s, base=2) for s in np.apply_along_axis(lambda x: ''.join(map(str, map(int, x))), 1, flags_check)] 
    dataRaw['flag_code'] = flag_code                                            # generate a binary representation of the flags
    print(dataRaw)                                                              # displays data so you can check its the timestep you wanted

    # resample clean data to 'StepSize' timestep
    cleanCheck = dataClean.copy()
    numeric = dataClean.select_dtypes(include=np.number).columns
    string = dataClean.select_dtypes(include=object).columns
    agg_dict = {**{c: 'mean' for c in numeric},**{c: 'first' for c in string}}
    dataClean = dataClean.resample(StepSize).agg(agg_dict)  # averages the numbers and returns the first string over the designated time step
    print('clean')
    num15 = dataClean[numeric].resample(StepSize).quantile(0.1587)
    txt15 = dataClean[string].resample(StepSize).first()
    data15 = num15.join(txt15)  
    print('15%')
    num84 = dataClean[numeric].resample(StepSize).quantile(0.8413)
    txt84 = dataClean[string].resample(StepSize).first()
    data84 = num84.join(txt84)  
    print('84%')
    mask = ((dataClean == 0) | (dataClean.isna())).all(axis=1)
    dataClean = dataClean.loc[~mask]
    data15 = data15.loc[~mask]
    data84 = data84.loc[~mask]
    times = dataClean.index.to_numpy()
    completeness = []
    for i in range(len(times)):
        ts = times[i]
        tf = ts +pd.Timedelta(f'1{StepSize}') - pd.Timedelta(1,'min')
        dt = pd.Timedelta(f'1{StepSize}').total_seconds()/60/3 #1 scan every 3 minutes
        slct = cleanCheck.loc[ts:tf]
        completeness.append(len(slct)/dt) #how many minutes/vs minutes in an hour
    dataClean['completeness'] = completeness

    dataSTP = ConvertToSTP(dataClean.copy())                                           # convert the hourly averaged cleaned data to stp
    data15 = ConvertToSTP(data15.copy())
    data84 = ConvertToSTP(data84.copy())
    num_cols= [s for s in data15.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
    data15 = data15[num_cols]
    data84 = data84[num_cols]
    data15 = data15.add_prefix('15.87% ', axis='columns')
    data84 = data84.add_prefix('84.13% ', axis='columns')
    data_EBAS = dataSTP.copy()
    data_EBAS['T_int'] = data_EBAS['Sheath Temp (C)'].to_numpy()+273.15
    data_EBAS['p_int'] = data_EBAS['Sheath Pressure (kPa)'].to_numpy()*10
    data_EBAS['RH_int'] = data_EBAS['Sheath Relative Humidity (%)'].to_numpy()
    kept_cols = ['completeness','T_int', 'p_int','RH_int']
    num_cols= [s for s in data_EBAS.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
    kept_cols.extend(num_cols)
    data_EBAS = data_EBAS[kept_cols]
    EBAS_cols = ['T_int', 'p_int','RH_int']
    for num in num_cols:
        n = data_EBAS.pop(num)
        data_EBAS[num] = n
        data_EBAS[f'15.87% {num}'] = data15[f'15.87% {num}'].to_numpy()
        data_EBAS[f'84.13% {num}'] = data84[f'84.13% {num}'].to_numpy()

    # Saves data upon request
    if dev and (out_name!=''): 
        raw_path = Path(inFolderpath).with_name(out_name+'.csv')
        dataRaw.to_csv(raw_path)
        clean_path = Path(inFolderpath).with_name(out_name+'_clean.csv')
        dataClean.to_csv(clean_path)
        stp_path = Path(inFolderpath).with_name(out_name+'_clean_stp.csv')
        dataSTP.to_csv(stp_path)
        print(f'3 files generated at {raw_path.parent}')
    elif not dev: 
        if (isinstance(Data, str))|(isinstance(Data, PurePath)):                # If there was a filepath passed, prompt for a name for the file, and save in the same folder as the filepath
            if folder == '': folder = Path(Data).parent
            name = input('\nEnter the desired name of your final files. [DO NOT INCLUDE .CSV]:\n' \
                        f'(This will place the file within {folder})\n')
            dataRaw.to_csv(Path(Data).with_name(name+'.csv')) 
            dataClean.to_csv(Path(Data).with_name(name+'_clean.csv')) 
            dataSTP.to_csv(Path(Data).with_name(name+'_clean_stp.csv')) 
            print(f'3 files generated at {folder}')
        else: 
            if folder =='':
                path_out = input('\nEnter the full path for your final files. [DO NOT INCLUDE .CSV]:\n')  
            else: 
                name = input('\nEnter the desired name of your final files. [DO NOT INCLUDE .CSV]:\n' \
                        f'(This will place the file within {folder})\n')
                path_out = folder+'/'+name
            dataRaw.to_csv(Path(path_out+'.csv')) 
            dataClean.to_csv(Path(path_out+'_clean.csv')) 
            dataSTP.to_csv(Path(path_out+'_clean_stp.csv')) 
            print(f'3 files generated at {path_out.parent}')
    return dataRaw, dataClean, dataSTP, data_EBAS                             # Return QA dataframe

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
    data : [Pandas DataFrame] outliers adding into original data frame
    """
    # generate a forward mean
    forward_mean = (data[name].shift(-1)        # Shift all rows back by 1
                    .rolling(window=size-1)     # Look forward at the next N-1 rows
                    .mean()                     # Mean these next N-1 rows (mean automatically placed at the index of the last row)
                    .shift(-(size-2)))          # Shift mean back from last index to first index
    # Check if the value of each scan exceeds the mean value of the next N-1 scans by more than a certain amount
    outliers = ((data[name] - forward_mean).abs() >avg_mult * forward_mean) 
    if name_out == '': name_out = name + ' flag'# Generate flag name if not provided
    data[name_out] = outliers.astype(int)
    return data   #return dataframe with the outlier added in as an additonal column

def FindOutliersRange(data,name,Min,Max,name_out = '',):
    """
    Takes in a data frame of TSI 3938 SMPS data and identifies outliers 
    within the data set that are outside of the Min and Max indicated values. 
    Returns a dataframe with the outlier flagged.
    ----------
    Parameters
    ++++++++++
    data : [Pandas DataFrame] SMPS dataframe with metadata removed
    name : [str] name of the column used to id outliers
    Max : [float] upper bounds for conditional flagging
    Min : float] lower bounds for conditional flagging
    name_out : [str] name of the flag column for outputting (default = '')

    Returns
    ++++++++++
    data : [Pandas DataFrame] outliers adding into original data frame
    """
    if name_out =='': name_out = name + ' flag'        # Generate flag name if not provided
    data[name_out] = ((data[name].to_numpy()<Min)|(data[name].to_numpy()>Max)).astype(int)
    return data  #return the dataframe with addition of flag

def ConvertToSTP(Data, Tstp = 273.15, Pstp= 101.325):
    if (isinstance(Data, str))|(isinstance(Data, PurePath)): # if data is a file path, read in Data to dataframe, otherwise assume passed dataframe
        Data = pd.read_csv(Data)
        Data = Data.set_index('Datetime(UTC)')
    cols = [col for col in Data.columns.to_numpy() if ('.' in col) and (col.split('.')[0].isdigit())] # pull out numeric columns to correct
    num_cols = cols.copy()
    cols.append('Total Concentration (#/cm³)') #stp correct Total Conc as well
    Tact = Data['Aerosol Temperature (C)'].to_numpy() +273.15
    Pact = Data['Sheath Pressure (kPa)'].to_numpy()
    Data['Standard Temperature (C)'] = Tstp - 273.15
    Data['Standard Pressure (kPa)'] = Pstp
    if conc =='num': # if numeric
        Tot = Data.pop('Total Concentration (#/cm³)')
        Data['Total Concentration (#/cm³)'] = Tot*Pstp/Pact*Tact/Tstp
    elif conc =='vol': # if volumetric
        Tot = Data.pop('Total Concentration (nm³/cm³)')
        Data['Total Concentration (nm³/cm³)'] = Tot*Pstp/Pact*Tact/Tstp
    elif conc =='mass': # if volumetric
        Tot = Data.pop('Total Concentration (µg/m³)')
        Data['Total Concentration (µg/m³)'] = Tot*Pstp/Pact*Tact/Tstp
    for col in num_cols: # move numeric cols to the end of Dataframe
        c = Data.pop(col)
        Data[col] = c*Pstp/Pact*Tact/Tstp
    return Data

def FinalCleaning(Data):
    if (isinstance(Data, str))|(isinstance(Data, PurePath)): # if data is a file path, read in Data to dataframe
        df = pd.read_csv(Data)
    elif isinstance(Data, pd.DataFrame): # If Data is a dataframe, reset index for processing
        df = Data.reset_index()
    lens = len(df.index)
    try:
        df = df.set_index('DateTime Sample Start')
        df.index = df.index.rename('Datetime(UTC)') # utilize a simpler index header common with the CCN 
        df.index = pd.to_datetime(df.index, format='mixed') 
    except:
        df = df.set_index('Datetime(UTC)')
        df.index = pd.to_datetime(df.index, format='mixed') 
    flags = [col for col in df.columns.to_numpy() if ('Flag' in col)]
    df = df.loc[(df[flags] == 0).all(axis=1)]                    # drop where any flag == 1
    drop = [col for col in df.columns.to_numpy() if ('raw' in col)|('_' in col)|('Flag' in col)]
    df = df.drop(columns=drop)                                   # drop the QA columns denoted with a _[diameter]
    nums = [s for s in df.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())&('raw' not in s)]
    slct = df[nums]
    mask = ((slct == 0) | (slct.isna())).all(axis=1)            # drop empty rows or all 0
    for date in bad_dates:                                      # remove periods of suspected bad data
        mask |= (df.index >= date[0]) & (df.index <= date[-1])
    df = df.loc[~mask]
    lenf = len(df.index.to_numpy())
    print(f'change = {lens}-{lenf}')
    return df

def SMPS_EBAS2(df,folder_out):
    """
    Takes in a path to processed SMPS file and generates a NASA AMES formated file
    ----------
    Parameters
    ++++++++++
    df : [panda Dataframe] processed SMPS file
    folder_out : [str/path-like] Path to folder to place EBAS file

    Returns
    ++++++++++
    NONE
    """
    # df=df.set_index('Datetime(UTC)') 
    df = df.fillna(0)
    df.index = pd.to_datetime(df.index)

    dates= df.index.to_list()
    # input(np.isnan(np.sum(df.values.tolist())))
    completeness = df['completeness'].values.tolist()
    num_cols = [s for s in df.columns.to_numpy() if ('%' not in s)&(('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
    num_cols = sorted(num_cols, key=lambda x: float(x))
    smps_header = ['pressure, hPa, Location=instrument internal sheath, Statistics=arithmetic mean, Matrix=instrument',
                   'relative_humidity, %, Location=instrument internal sheath, Statistics=arithmetic mean, Matrix=instrument',
                   'temperature, K, Location=instrument internal, Statistics=arithmetic mean, Matrix=instrument']
    smps_units = ['K','hPa', '%']
    smps_stats = ['Statistics=arithmetic mean','Statistics=arithmetic mean','Statistics=arithmetic mean']
    for num in num_cols:
        smps_header.append(f'particle_number_size_distribution, 1/cm3, D={num} nm, Statistics=arithmetic mean')
        smps_header.append(f'particle_number_size_distribution, 1/cm3, D={num} nm, Statistics=percentile:15.87')
        smps_header.append(f'particle_number_size_distribution, 1/cm3, D={num} nm, Statistics=percentile:84.13')
    smps_cols = ['T_int', 'p_int','RH_int']
    data_cols = smps_cols.copy()
    for num in num_cols:
        data_cols.append(num)
        data_cols.append(f'15.87% {num}')
        data_cols.append(f'84.13% {num}')
    smps_concs = []
    for num in num_cols:
        smps_concs.append(f'Conc {num}')
        smps_concs.append(f'15prc {num}')
        smps_concs.append(f'84prc {num}')
    smps_cols.extend(smps_concs)
    df = df.drop(columns = ['completeness'])
    data = df[data_cols].values.tolist()

    # df['flag'] = [000]
    # df['flag']['Q_flag' ==1] = [662]
    # df['flag']['integrity_flag' == 1] = [111]
    flags = [[[000] for j in range(len(data[i]))] for i in range(len(dates))]
    for i in range(len(completeness)):
        if (completeness[i] < .9)&(completeness[i]>0.75):
            flags[i] = [[394] for j in range(len(data[i]))]
        elif (completeness[i] < .75)&(completeness[i]>0.5):
            flags[i] = [[392] for j in range(len(data[i]))]
        elif (completeness[i] < .5):
            flags[i] = [[390] for j in range(len(data[i]))]
    data =  list(map(list, zip(*data)))
    flags = list(map(list, zip(*flags)))

    date_list = []
    for i in range(len(dates)-1):
        if i == 0:
            date_list.append(pd.to_datetime(dates[i]))
        elif (pd.to_datetime(dates[i+1])-pd.to_datetime(dates[i]) > pd.Timedelta(1, 'hr')):
            date_list.append(pd.to_datetime(dates[i]))
            date_list.append(pd.to_datetime(dates[i+1]))
    date_list.append(pd.to_datetime(dates[-1]))  
    # data = np.nan_to_num(data, nan=1, posinf=1e6, neginf=-1e6)
    # input(np.isnan(np.sum(data)))
    ebas_genfile(folder_out, data, flags, dates, smps_header, smps_cols)

def SMPS_EBAS(df,folder_out):
    """
    Takes in a path to processed SMPS file and generates a NASA AMES formated file
    ----------
    Parameters
    ++++++++++
    df : [panda Dataframe] processed SMPS file
    folder_out : [str/path-like] Path to folder to place EBAS file

    Returns
    ++++++++++
    NONE
    """
    # df=df.set_index('Datetime(UTC)') 
    df = df.fillna(0)
    df.index = pd.to_datetime(df.index)

    dates= df.index.to_list()
    # input(np.isnan(np.sum(df.values.tolist())))
    completeness = df['completeness'].values.tolist()
    num_cols = [s for s in df.columns.to_numpy() if ('%' not in s)&(('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
    num_cols = sorted(num_cols, key=lambda x: float(x))
    smps_header = ['temperature','pressure','relative_humidity',]
    smps_units = ['K','hPa', '%']
    smps_stats = ['arithmetic mean','arithmetic mean','arithmetic mean']
    smps_locs = ['instrument internal sheath','instrument internal sheath','instrument internal sheath']
    smps_D = ['','','']
    for num in num_cols:
        smps_header.append(f'particle_number_size_distribution')
        smps_units.append('1/cm3')
        smps_D.append(str(num))
        smps_stats.append('arithmetic mean')
        smps_locs.append('')
        smps_header.append(f'particle_number_size_distribution')
        smps_units.append('1/cm3')
        smps_D.append(str(num))
        smps_stats.append('percentile:15.87')
        smps_locs.append('')
        smps_header.append(f'particle_number_size_distribution')
        smps_units.append('1/cm3')
        smps_D.append(str(num))
        smps_stats.append('percentile:84.13')
        smps_locs.append('')
    smps_cols = ['T_int', 'p_int','RH_int']
    data_cols = smps_cols.copy()
    for num in num_cols:
        data_cols.append(num)
        data_cols.append(f'15.87% {num}')
        data_cols.append(f'84.13% {num}')
    smps_concs = []
    for num in num_cols:
        smps_concs.append(f'Conc {num}')
        smps_concs.append(f'15prc {num}')
        smps_concs.append(f'84prc {num}')
    smps_cols.extend(smps_concs)
    df = df.drop(columns = ['completeness'])
    data = df[data_cols].values.tolist()

    # df['flag'] = [000]
    # df['flag']['Q_flag' ==1] = [662]
    # df['flag']['integrity_flag' == 1] = [111]
    flags = [[[000] for j in range(len(data[i]))] for i in range(len(dates))]
    for i in range(len(completeness)):
        if (completeness[i] < .9)&(completeness[i]>0.75):
            flags[i] = [[394] for j in range(len(data[i]))]
        elif (completeness[i] < .75)&(completeness[i]>0.5):
            flags[i] = [[392] for j in range(len(data[i]))]
        elif (completeness[i] < .5):
            flags[i] = [[390] for j in range(len(data[i]))]
    data =  list(map(list, zip(*data)))
    flags = list(map(list, zip(*flags)))
    ebas_genfile(folder_out, data, flags, dates, smps_header, smps_units, smps_stats, smps_locs, smps_D, smps_cols)

#if the program exists, run it
if __name__:
    main()