'''
Ben Sykes
4/19/2026
Description:
Reindexed AQS master files generated using the download_aqs.py and combine_distance.py processing. Pass
file into the datetime_adjust() function and the output files are automatically processed and outputted 
to the same directory as the inputted file. 
'''
import pandas as pd
import numpy as np
from os.path import basename, join, dirname

def datetime_adjust(f):
    """
    Takes in a path to yearly AQS file output. Reindex for datetime and coords. Output a continuous
    PM2.5 CSV, a daily PM2.5 CSV and a daily speciated CSV to the file path of the original AQS file.
    ----------
    Paramaters
    ++++++++++
    f : [str/path-like] Path to the AQS file

    Returns
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    """
    #Read in file, split year from original file and dir for output files 
    data= pd.read_csv(f)
    dir_name = dirname(f)
    file_name = basename(f)
    yearcsv = list(file_name.split('_'))[-1]
    year = yearcsv.replace('.csv', '')
    #Combine date and time into a single datetime column
    data['Datetime(UTC)'] = pd.to_datetime(data['date_gmt'] + ' ' + data['time_gmt'])
    #Combine lat and long coords into a single possition column
    data['Position'] = np.array(list(map(str, data["latitude"].to_numpy()))) +',' + np.array(list(map(str, data["longitude"].to_numpy())))
    #Combine county and state into single location column
    data['Location'] = np.array(list(map(str,data['county'])))+' county,'+np.array(list(map(str,data['state'])))
    #Unit conversion dictionary to increase readability
    #[LC stands for Local Conditions(ATP)]
    units = {
        "Micrograms/cubic meter (LC)": '[ug/m3 ATP]', 
        "Micrograms/cubic meter (25 C)": '[ug/m3 STP]',
        "Inverse 100 Megameters" : '[1/100*Mm]', 
        "Parts per billion": '[ppb]',
        "Parts per million": '[ppm]',
        'Millimeters (mercury)':'[mmHg]',
        'Degrees Centigrade':'[C]'}
    # Combine paramter and updated units into single column header
    label = [f"{a} {b}" for a, b in zip(data['parameter'].to_numpy(),np.array([units[key] for key in data['units_of_measure'].to_numpy() if key in units]))]
    data['label'] = label
    # Clean up data to only use needed columns
    # Use date and position as index variables
    data = data.set_index(['Datetime(UTC)','Position'])
    kept = ['Location','label','sample_measurement','uncertainty','qualifier','sample_frequency','method', 'dist_km'] #observation_percent
    bad_cols = [col for col in data.columns.to_numpy() if col not in kept]
    data = data.drop(columns =bad_cols)
    data = data[kept]
    data = data.sort_index()
    #Filter out lower frequency datasets from continuous 
    total = data[(data['sample_frequency']!='HOURLY')]
    cont = data[data['sample_frequency']=='HOURLY']
    #clean and pivot 3 seperatre datasets
    total = total.drop(columns='sample_frequency')
    cont = cont.drop(columns='sample_frequency')
    #Process Continuous datasets
    cont = cont.reset_index()
    cont = pd.pivot_table(cont, values='sample_measurement', index=['Datetime(UTC)'],
                       columns=['Position'], aggfunc="mean")
    cont['Measurement'] ='PM2.5 [ug/m3 ATP]' # all continous measurements are of PM2.5
    #Processed non-hourly Re-indexed 
    total = total.reset_index()
    total['Date(UTC)'] = total['Datetime(UTC)'].dt.date
    total = pd.pivot_table(total, values='sample_measurement', index=['Date(UTC)','Position','Location',],
                       columns=['label'], aggfunc="mean")
    # Output CSVs of updated data
    # Clean up columns for outputted file
    total.columns = (total.columns
                    .str.replace('TOR ', '', regex=False)
                    .str.replace('LC ', '', regex=False)
                    .str.replace('Ammonium Ion','NH4', regex=False)
                    .str.replace('Sulfate','SO4', regex=False)
                    .str.replace('Total Nitrate','NO3', regex=False)
                    .str.replace(' - Local Conditions', '', regex=False)
                    .str.replace('Average Ambient Pressure ', 'P_act', regex=False)
                    .str.replace('Average Ambient Temperature ', 'T_act', regex=False))
    # Apply STP conv
    total = STP_conv(total)
    # Calculated mass fractions for speciated data
    total = Frac_conv(total)
    # Drop negative values
    num = total._get_numeric_data()
    num[num<0] = 0
    # Split total dataframe into a speciated dataset and a PM2.5 dataset
    tot_cols= total.columns.to_numpy()
    for col in [c for c in tot_cols if 'OC' in c]:
        new_col = col.replace('OC', 'Org')
        pos = total.columns.get_loc(col)
        total.insert(pos + 1, new_col, total[col].to_numpy()*2) #apply conversion from OC to Org Matter 
    tot_cols= total.columns.to_numpy()
    chems= ['OC', 'Org', 'EC', 'NO3', 'NH4', 'SO4']
    conv = ['P_act','P_std','T_act','T_std']
    chem_cols = [col for col in tot_cols if (any(y in col for y in chems))|(any(y in col for y in conv))] #<-- Speciated columns 
    PM25_cols = [col for col in tot_cols if (not any(y in col for y in chems))&('PM2.5' in col)|(any(y in col for y in conv))] #<-- PM2.5 columns
    # Split dataframes and drop empty rows
    spec = total[chem_cols]
    PM25 = total[PM25_cols]
    PM25 = PM25.dropna(how ='all')
    spec = spec.dropna(thresh=5)
    # Output datafiles
    total.to_csv(join(dir_name,f"AQS_Reindexed_{year}.csv"))
    spec.to_csv(join(dir_name,f"AQS_combined_speciated_{year}.csv"))
    PM25.to_csv(join(dir_name,f"AQS_combined_PM25_{year}.csv"))
    cont.to_csv(join(dir_name,f"AQS_combined_PM25_cont_{year}.csv"))
    print('done')
    return data,spec

def STP_conv(data, sT = 0, sP= 760):
    """
    Takes in a AQS dataframe. Applies a STP conversion
    Cstp = Pstp/Pact * Tact/Tstp * Cact
    ----------
    Paramaters
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    sT : [float] standard temperature [C] (default = 0)
    sP : [float] standard pressure [mmHg] (default = 760)

    Returns
    ++++++++++
    data : [DataFrame] AQS file with STP conversions
    """
    cols = list(data.columns.values)            # List of all of the columns in the df
    cols.pop(cols.index('T_act[C]')) 
    cols.pop(cols.index('P_act[mmHg]'))
    data = data[cols+['T_act[C]','P_act[mmHg]']] # Move T and P values to end of dataframe
    # Add STP conversion values to end of dataframe
    data['P_std[mmHg]'] = sP
    data['T_std[C]'] = sT
    # DisplaY STP conv value for concentrations such that C_STP = STP_Conv * C_ATP 
    data['STP_Conv'] = sP/data['P_act[mmHg]'] * (data['T_act[C]']+273.15)/(sT+273.15)
    cols= [c for c in data.columns.to_numpy() if "ug/m3" in c]
    for col in cols:
        new_col = col.replace("ATP", "STP")
        data[new_col] = data[col]*data['STP_Conv']
    return data

def Frac_conv(data, pm_col = 'PM2.5 [ug/m3 ATP]',chems= ['OC', 'EC', 'NO3', 'NH4', 'SO4']):
    """
    Takes in a AQS dataframe. Converts mass concentrations to mass fractions.
    ----------
    Paramaters
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    pm_col : [str] PM2.5 column for total mass
    chems : [list of str] Species to calculate mass fraction for

    Returns
    ++++++++++
    data : [DataFrame] AQS file with mass fraction conversions
    """
    cols= [c for c in data.columns.to_numpy() if (any(y in c for y in chems))&('ATP' in c)]
    for col in cols:
        new_col = col.replace(' PM2.5 [ug/m3 ATP]', "/total")
        data[new_col] = data[col]/data[pm_col]
    return data

if __name__ == "__main__":
    f = r"C:\Users\bensy\Documents\Research\AQS_CSVs_for_Reindexing\AQS_combined_within_350km_2025.csv" #input('Input path to yearly combined AQS file: ')
    d,s= datetime_adjust(f)