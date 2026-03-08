import pandas as pd
import numpy as np
from os.path import basename, join, dirname

def datetime_adjust(f):
    """
    Takes in a path to yearly AQS file output. Reindex for datetime and coords. Output a continuous
    PM2.5 CSV, a daily PM2.5 CSV and a daily speciated CSV to the file path of the original AQS file
    ----------
    Paramaters
    ++++++++++
    path : [str/path-like] Path to the AQS file

    Returns
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    """
    #Read in file, split year from original file and dir to output files too
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
    #Unit conversion dictionary to increase readability/ LC stands for Local Conditions(ATP), but local conditions not dir specified
    units = {
        "Micrograms/cubic meter (LC)": '[µg/m^3 ATP]',
        "Micrograms/cubic meter (25 C)": '[µg/m^3 STP]',
        "Inverse 100 Megameters" : '[1/100*Mm]',
        "Parts per billion": '[ppb]',
        "Parts per million": '[ppm]'}
    #Combine paramter and updated units into single column header
    label = [f"{a} {b}" for a, b in zip(data['parameter'].to_numpy(),np.array([units[key] for key in data['units_of_measure'].to_numpy() if key in units]))]
    data['label'] = label
    #clean up data to only use needed columns
    data = data.set_index(['Datetime(UTC)','Position'])
    kept = ['Location','label','sample_measurement','uncertainty','qualifier','sample_frequency','method', 'dist_km'] #observation_percent
    cols = [col for col in data.columns.to_numpy() if col not in kept]
    data = data.drop(columns =cols)
    data = data[kept]
    data = data.sort_index()
    #Filter out speciated and two PM2.5 datasets
    spec = data[(data['sample_frequency']!='HOURLY')&(data['label']!='PM2.5 - Local Conditions [µg/m^3 ATP]')]
    PM25 = data[(data['sample_frequency']!='HOURLY')&(data['label']=='PM2.5 - Local Conditions [µg/m^3 ATP]')]
    cont = data[data['sample_frequency']=='HOURLY']
    #clean and pivot 3 seperatre datasets
    spec = spec.drop(columns='sample_frequency')
    PM25 = PM25.drop(columns='sample_frequency')
    cont = cont.drop(columns='sample_frequency')
    #Cont
    cont = cont.reset_index()
    cont = pd.pivot_table(cont, values='sample_measurement', index=['Datetime(UTC)'],
                       columns=['Position'], aggfunc="mean")
    cont['Measurement'] ='PM2.5 [µg/m^3 ATP]'
    #Daily PM2.5
    PM25 = PM25.reset_index()
    PM25['Date(UTC)'] = PM25['Datetime(UTC)'].dt.date
    PM25 = pd.pivot_table(PM25, values='sample_measurement', index=['Date(UTC)','Position','Location',],
                       columns=['label'], aggfunc="mean")
    #Speciated
    spec = spec.reset_index()
    spec['Date(UTC)'] = spec['Datetime(UTC)'].dt.date
    spec = pd.pivot_table(spec, values='sample_measurement', index=['Date(UTC)','Position','Location',],
                       columns=['label'], aggfunc="mean")
    #Output CSVs of updated data
    data.to_csv(join(dir_name,f"AQS_combined_{year}.csv"))
    spec.to_csv(join(dir_name,f"AQS_combined_speciated_{year}.csv"))
    PM25.to_csv(join(dir_name,f"AQS_combined_PM25_{year}.csv"))
    cont.to_csv(join(dir_name,f"AQS_combined_cont_{year}.csv"))
    return data

if __name__ == "__main__":
    f = input('Input path to yearly combined AQS file: ')
    datetime_adjust(f)
