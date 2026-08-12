'''
Ben Sykes
5/11/2026
Description:
Reindexed AQS master files generated using the download_aqs.py and combine_distance.py processing. Pass
file into the datetime_adjust() function and the output files are automatically processed and outputted 
to the same directory as the inputted file. 
'''
import pandas as pd
import numpy as np
from os.path import basename, join, dirname
import os

def county_conv(cnty):
    '''Provided a county name -> return whether that county is rural or Urban
    ----------
    Paramaters
    ++++++++++
    cnty : [str] county name as key for dictionary 

    Returns
    ++++++++++
    *val : [str] dictionary value associated with the passed key {"Urban" or "Rural"}
    '''
    conv_dict = {"Alamance county,North Carolina": "Urban",
            "Albemarle county,Virginia": "Urban",
            "Amherst county,Virginia": "Rural",
            "Avery county,North Carolina": "Rural",
            "Bell county,Kentucky": "Rural",
            "Blount county,Tennessee": "Urban",
            "Boyd county,Kentucky": "Rural",
            "Bristol City county,Virginia": "Rural",
            "Buchanan county,Virginia": "Rural",
            "Buncombe county,North Carolina": "Urban",
            "Cabell county,West Virginia": "Urban",
            "Carter county,Kentucky": "Rural",
            "Caswell county,North Carolina": "Rural",
            "Catawba county,North Carolina": "Urban",
            "Chatham county,North Carolina": "Urban",
            "Chesterfield county,South Carolina": "Rural",
            "Clarke county,Georgia": "Rural",
            "Cumberland county,North Carolina": "Urban",
            "Davidson county,North Carolina": "Rural",
            "Durham county,North Carolina": "Urban",
            "Edgefield county,South Carolina": "Rural",
            "Fayette county,Kentucky": "Urban",
            "Florence county,South Carolina": "Urban",
            "Forsyth county,North Carolina": "Urban",
            "Gaston county,North Carolina": "Urban",
            "Greenville county,South Carolina": "Urban",
            "Guilford county,North Carolina": "Urban",
            "Gwinnett county,Georgia": "Urban",
            "Hall county,Georgia": "Urban",
            "Hamilton county,Tennessee": "Urban",
            "Haywood county,North Carolina": "Rural",
            "Jackson county,North Carolina": "Rural",
            "Jessamine county,Kentucky": "Urban",
            "Johnston county,North Carolina": "Urban",
            "Kanawha county,West Virginia": "Urban",
            "Knox county,Tennessee": "Urban",
            "Lee county,North Carolina": "Urban",
            "Lexington county,South Carolina": "Urban",
            "Loudon county,Tennessee": "Rural",
            "Lynchburg City county,Virginia": "Rural",
            "Madison county,Kentucky": "Rural",
            "McDowell county,North Carolina": "Rural",
            "McMinn county,Tennessee": "Rural",
            "Mecklenburg county,North Carolina": "Urban",
            "Mitchell county,North Carolina": "Rural",
            "Montgomery county,North Carolina": "Rural",
            "Montgomery county,Virginia": "Rural",
            "Murray county,Georgia": "Rural",
            "Oconee county,South Carolina": "Rural",
            "Perry county,Kentucky": "Rural",
            "Pike county,Kentucky": "Rural",
            "Pulaski county,Kentucky": "Rural",
            "Putnam county,Tennessee": "Rural",
            "Raleigh county,West Virginia": "Urban",
            "Richland county,South Carolina": "Urban",
            "Richmond county,Georgia": "Urban",
            "Roane county,Tennessee": "Rural",
            "Roanoke City county,Virginia": "Urban",
            "Roanoke county,Virginia": "Rural",
            "Robeson county,North Carolina": "Rural",
            "Rockbridge county,Virginia": "Rural",
            "Rowan county,North Carolina": "Rural",
            "Russell county,Kentucky": "Rural",
            "Salem City county,Virginia": "Rural",
            "Spartanburg county,South Carolina": "Urban",
            "Sullivan county,Tennessee": "Urban",
            "Swain county,North Carolina": "Rural",
            "Wake county,North Carolina": "Urban",
            "Watauga county,North Carolina": "Rural",
            "Wayne county,North Carolina": "Urban",
            "Wood county,West Virginia": "Rural",
            "York county,South Carolina": "Urban"}
    if cnty in list(conv_dict.keys()):
        return conv_dict[cnty]
    else:
        return 'Void'

def datetime_adjust(chem, pm25, out, freq ='ME', freq_total = 'W', state_county = 0):
    """
    Takes in a path to yearly AQS file output. Reindex for datetime and coords. Output a continuous
    PM2.5 CSV, a daily PM2.5 CSV and a daily speciated CSV to the file path of the original AQS file.
    ----------
    Paramaters
    ++++++++++
    chem : [str/path-like] Path to the AQS chemistry file
    pm25 : [str/path-like] Path to the AQS pm2.5 file
    out : [str/path-like] Output Path for the combined file
    freq : [str] Frequency to process combined datafile to (default = ME)
    freq_total : [str] Frequency to process reindex datafile to (default = W)
    state_county : [list-like/0] Only reports data from the selected state and county (default = 0, not used)

    Returns
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    spec : [DataFrame] Speciated chemistry data from the AQS re-indexed file
    Combined : [DataFrame] Rural|Urban|Total combined AQS dataset 
    nearby : [DataFrame] Nearby(avery county) AQS dataset
    """
    #Read in file, split year from original file and dir for output files 
    no_chem = False
    no_mass = False
    print(chem)
    if chem == '': no_chem = True
    if pm25 =='': no_mass = True
    if not no_chem:
        chem_data= pd.read_csv(chem)
    if not no_mass:
        pm25_data= pd.read_csv(pm25)
    dir_name = out
    file_name = basename(pm25)
    if no_chem:
        data = pm25_data
    elif no_mass:
        data = chem_data
    else:
        data = pd.concat([chem_data, pm25_data], ignore_index=True)
    if state_county != 0:
        data = data[data['state_code'] == (state_county[0])]
        data = data[data['county_code'] == (state_county[-1])]
    # if '2026' in pm25:
    #     input(data)
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
    kept = ['Location','Category','label','sample_measurement','uncertainty','qualifier','sample_frequency','method', 'dist_km'] #observation_percent
    bad_cols = [col for col in data.columns.to_numpy() if col not in kept]
    data['Category'] = [county_conv(loc) for loc in data['Location'].to_numpy()]
    data = data.drop(columns =bad_cols)
    data = data[kept]
    data = data.sort_index()
    #Filter out lower frequency datasets from continuous 
    total = data[(data['sample_frequency']!='HOURLY')]
    cont = data[data['sample_frequency']=='HOURLY']
    #clean and pivot 2 seperatre datasets
    total = total.drop(columns='sample_frequency')
    cont = cont.drop(columns='sample_frequency')
    #Process Continuous datasets to generate a continuous PM2.5 measurement for filling in gaps of PM2.5 data
    cont = cont.reset_index()
    if not cont.empty:
        cont['columns'] = [f"{a}|{b}" for a, b in zip(cont['Location'].to_numpy(),cont['Position'].to_numpy())]
        cont = pd.pivot_table(cont, values='sample_measurement', index=['Datetime(UTC)'],
                        columns=['columns'], aggfunc="mean")
        cont.index =  pd.to_datetime(cont.index, format='mixed')
        cont = cont.clip(lower=0)
        try:
            cont = cont.drop(columns=['Blount county,Tennessee|35.63348,-83.941606'])
        except:
            print('Blount county,Tennessee cont. PM measurment not found')
        daily = cont.resample('d').mean()
        cont['Measurement'] ='PM2.5 [ug/m3 ATP]' # all continous measurements are of PM2.5

        daily_long = daily.stack().reset_index()
        daily_long['Datetime(UTC)'] = daily_long['Datetime(UTC)'].dt.date
        daily_long.columns = ['Date(UTC)', "poc", "PM2.5 cont.[ug/m3 ATP]"]
        loc_split = daily_long["poc"].str.split("|", expand=True)

        daily_long["Location"] = loc_split[0]
        daily_long["Position"] = loc_split[1]

        daily_long = daily_long.drop(columns=["poc"])
    #Processed non-hourly Re-indexed 
    total = total.reset_index()
    total['Date(UTC)'] = pd.to_datetime(total['Datetime(UTC)'].dt.date)
    total = pd.pivot_table(total, values='sample_measurement', index=['Date(UTC)','Position','Category','Location'],
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
                    .str.replace('Average Ambient Temperature ', 'T_act', regex=False)
                    .str.replace('Acceptable PM2.5 AQI & Speciation Mass ', 'PM2.5 accept', regex=False))
    # Add back cont. pm2.5 for processing
    if not cont.empty:
        total = Cont_combo(total, daily_long)
    else: 
        total = total.reset_index()
        total = total.set_index(["Date(UTC)", "Location", "Position"])
    #Resample to chosen resolution, usually weekly due to 3-6 Day sample frequency 
    if freq_total != 'D':
        try:
            total['T_act[C]'] = (total.groupby(level=['Location', 'Position'])['T_act[C]'].transform(lambda s: s.interpolate(method='linear')))
            total['P_act[mmHg]'] = (total.groupby(level=['Location', 'Position'])['P_act[mmHg]'].transform(lambda s: s.interpolate(method='linear')))
            total = total.groupby(level=['Location', 'Position']).ffill(limit=4)
        except:
            total = total.groupby(level=['Location', 'Position']).ffill(limit=4)
        total = (total.groupby(level=['Location', 'Position']).resample(freq_total, level='Date(UTC)')
                .agg({'Category': 'first'} | {col: 'mean'for col in total.select_dtypes('number').columns}))
        total = (total.reorder_levels(['Date(UTC)', 'Location', 'Position']).sort_index(level='Date(UTC)'))
    ''' Apply STP conversions '''
    try:
        total = STP_conv(total)
    except Exception as e: 
        input(f'Failed to preform StP conversions for {year} due to {e}')
    ''' Calculated mass fractions for speciated data '''
    if not no_chem:
        try:
            total = Frac_conv(total)
        except Exception as e: 
            input(f'Failed to preform Mass Fraction conversions for {year} due to {e}')
        #apply conversion from OC to Org Matter and insert column directly after OC
        tot_cols= total.columns.to_numpy()
        for col in [c for c in tot_cols if 'OC' in c]:
            new_col = col.replace('OC', 'Org')
            pos = total.columns.get_loc(col)
            total.insert(pos + 1, new_col, total[col].to_numpy()*2) 
    
    ''' Group and average using Rural vs Urban Distinctions along with a combined Total Column '''
    # Resample to a monthly frequency
    reset = total.reset_index()
    grouped = reset.groupby(["Date(UTC)", "Category"]).mean(numeric_only=True).reset_index()
    Urban = grouped[grouped['Category'] =='Urban']
    Rural = grouped[grouped['Category'] =='Rural']
    #Drop measurements reported at actual temperature and pressure
    drop = [col for col in Rural.columns.to_numpy() if ("ATP" in col) or ('std' in col) or ('act' in col) or ('Conv' in col)
             or ('cont.' in col) or ('accept' in col) or ('filter' in col)]
    Urban = Urban.drop(columns=['Category'])
    Urban = Urban.drop(columns=drop)
    Rural = Rural.drop(columns=['Category'])
    Rural = Rural.drop(columns=drop)
    Urban = Urban.set_index('Date(UTC)')
    Rural = Rural.set_index('Date(UTC)')
    Rural = Rural.add_prefix('Rural ')
    Urban = Urban.add_prefix('Urban ')
    cols_filt = []
    for col in Urban.columns.to_numpy():
        cols_filt.append(col)
        rural = col.replace("Urban", "Rural")
        cols_filt.append(rural)
    Combined = pd.merge(Urban,Rural, left_index=True, right_index=True)
    Combined = Combined[cols_filt]
    Combined.index = pd.to_datetime(Combined.index, format='mixed')
    Combined = Combined.resample(freq).mean()
    Combined.replace([np.inf, -np.inf], np.nan, inplace=True)
    ordered = []
    for urban in Urban.columns.to_numpy():
        ordered.append(urban)
        rural = urban.replace("Urban", "Rural")
        ordered.append(rural)
        tot = urban.replace("Urban", "Total")
        Combined[tot] = Combined[[urban, rural]].mean(axis=1)
        ordered.append(tot)
    Combined = Combined[ordered]

    ''' Pull out just the Avery Dataset '''
    # Resample to a monthly frequency
    reset = total.reset_index()
    reset = reset[reset['Location'] =='Avery county,North Carolina']
    nearby = reset.groupby(["Date(UTC)"]).mean(numeric_only=True)
    nearby.index = pd.to_datetime(nearby.index, format='mixed')
    drop = [col for col in nearby.columns.to_numpy() if ('cont.' in col) or ('accept' in col) or ('filter' in col)
            or ('P_std' in col) or ('T_std' in col)]
    nearby = nearby.drop(columns=drop)
    nearby = nearby.resample(freq).mean()
    nearby = nearby.dropna(axis=1, how='all')
    nearby.to_csv(join(dir_name,f"AQS_avery_{year}.csv"))
    # Drop negative values within the datasets
    num = total._get_numeric_data()
    num[num<0] = 0
    # Split total dataframe into a speciated dataset and a PM2.5 dataset
    tot_cols= total.columns.to_numpy()
    chems= ['OC', 'Org', 'EC', 'NO3', 'NH4', 'SO4']
    conv = ['P_act','P_std','T_act','T_std']
    #split into PM2.5 and speciated columns. Include conversion factors in both datasets, for future efforts
    chem_cols = [col for col in tot_cols if (any(y in col for y in chems))|(any(y in col for y in conv))] #<-- Speciated columns 
    PM25_cols = [col for col in tot_cols if (not any(y in col for y in chems))&('PM2.5' in col)|(any(y in col for y in conv))] #<-- PM2.5 columns
    # Split dataframes and drop empty rows
    spec = total[chem_cols]
    PM25 = total[PM25_cols]
    PM25 = PM25.dropna(how ='all')
    spec = spec.dropna(thresh=5)
    # Output datafiles
    freq_dict = {'ME': 'Monthly', 'W': 'Weekly', "d": 'Daily'}
    Combined.to_csv(join(dir_name,f"AQS_{freq_dict[freq]}_{year}.csv"))
    total.to_csv(join(dir_name,f"AQS_Reindexed_{year}.csv"))
    spec.to_csv(join(dir_name,f"AQS_combined_speciated_{year}.csv"))
    PM25.to_csv(join(dir_name,f"AQS_combined_PM25_{year}.csv"))
    cont.to_csv(join(dir_name,f"AQS_combined_PM25_cont_{year}.csv"))
    print(f'Finished for {year}')
    return data,spec,Combined, nearby

def large_adjust(large, out, freq ='ME', freq_total = 'd'):
    """
    Takes in a path to pre-processed yearly AQS file input. Reindex for datetime and coords. Output a
    daily CSV to the specified file path.
    ----------
    Paramaters
    ++++++++++
    chem : [str/path-like] Path to the AQS chemistry file
    pm25 : [str/path-like] Path to the AQS pm2.5 file
    out : [str/path-like] Output Path for the combined file
    freq : [str] Frequency to process combined datafile to (default = ME)
    freq_total : [str] Frequency to process reindex datafile to (default = W)
    state_county : [list-like/0] Only reports data from the selected state and county (default = 0, not used)

    Returns
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    spec : [DataFrame] Speciated chemistry data from the AQS re-indexed file
    Combined : [DataFrame] Rural|Urban|Total combined AQS dataset 
    nearby : [DataFrame] Nearby(avery county) AQS dataset
    """
    #Read in file, split year from original file and dir for output files 
    data= pd.read_csv(large)
    good_states = [37,13,21,45,47,51,54]
    data = data.loc[data['State Code'].isin(good_states)]
    dir_name = out
    file_name = basename(large)
    yearcsv = list(file_name.split('_'))[-1]
    year = yearcsv.replace('.csv', '')
    #Generate date column as a datetime
    try:
        data['Date(UTC)'] = pd.to_datetime(data['date_gmt'])
    except:
        data['Date(EST)'] = pd.to_datetime(data['Date Local'])
    #Combine lat and long coords into a single possition column
    data['Position'] = np.array(list(map(str, data["Latitude"].to_numpy()))) + ',' + np.array(list(map(str, data["Longitude"].to_numpy())))
    #Combine county and state into single location column
    data['Location'] = np.array(list(map(str,data['County Name'])))+' county,'+np.array(list(map(str,data['State Name'])))
    #Unit conversion dictionary to increase readability
    #[LC stands for Local Conditions(ATP)]
    units = {
        "Micrograms/cubic meter (LC)": '[ug/m3 ATP]', 
        "Micrograms/cubic meter (25 C)": '[ug/m3 STP]',
        "Inverse 100 Megameters" : '[1/100*Mm]', 
        "Parts per billion": '[ppb]',
        "Parts per billion Carbon" : '[ppbC]',
        "Parts per million": '[ppm]',
        'Millimeters (mercury)':'[mmHg]',
        'Degrees Centigrade':'[C]'}
    # Combine paramter and updated units into single column header
    label = [f"{a} {b}" for a, b in zip(data['Parameter Name'].to_numpy(),np.array([units[key] for key in data['Units of Measure'].to_numpy() if key in units]))]
    data['label'] = label
    # Clean up data to only use needed columns
    # Use date and position as index variables
    data = data.set_index(['Date(EST)','Position'])
    kept = ['Location','Category','label','Arithmetic Mean','1st Max Value'] #observation_percent
    bad_cols = [col for col in data.columns.to_numpy() if col not in kept]
    data['Category'] = [county_conv(loc) for loc in data['Location'].to_numpy()]
    data.loc[data['Category'] == 'Void'] = np.nan
    data = data.drop(columns =bad_cols)
    data = data[kept]
    data = data.sort_index()
    #Processed Daily Re-indexed 
    total = data.reset_index()
    total = pd.pivot_table(total, values='Arithmetic Mean', index=['Date(EST)','Position','Category','Location'],
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
                    .str.replace('Average Ambient Temperature ', 'T_act', regex=False)
                    .str.replace('Acceptable PM2.5 AQI & Speciation Mass ', 'PM2.5 accept', regex=False))
    total = total.reset_index()
    total = total.set_index(["Date(EST)", "Location", "Position"])
    #Resample to chosen resolution, usually weekly due to 3-6 Day sample frequency 
    if freq_total != 'D':
        total = total.groupby(level=['Location', 'Position']).ffill(limit=4)
        total = (total.groupby(level=['Location', 'Position']).resample(freq_total, level='Date(EST)')
                .agg({'Category': 'first'} | {col: 'mean'for col in total.select_dtypes('number').columns}))
        total = (total.reorder_levels(['Date(EST)', 'Location', 'Position']).sort_index(level='Date(EST)'))
    # Drop negative values within the datasets
    num = total._get_numeric_data()
    num[num<0] = 0
    total = total.dropna(axis='index', thresh=10)
    # Output datafiles
    print(dir_name)
    total.to_csv(join(dir_name,f"Reindexed_{file_name}.csv"))
    print(f'Finished for {year}')
    return data,total

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
    cols = list(data.columns.values)                   # List of all of the columns in the df
    cols.pop(cols.index('T_act[C]')) 
    cols.pop(cols.index('P_act[mmHg]'))
    data = data[cols+['T_act[C]','P_act[mmHg]']].copy() # Move T and P values to end of dataframe
    # Add STP conversion values to end of dataframe
    data['P_std[mmHg]'] = sP
    data['T_std[C]'] = sT
    # DisplaY STP conv value for concentrations such that C_STP = STP_Conv * C_ATP 
    data['STP_Conv'] = sP/data['P_act[mmHg]'] * (data['T_act[C]']+273.15)/(sT+273.15)
    cols= [c for c in data.columns.to_numpy() if "ug/m3" in c]
    new_data = {col.replace("ATP", "STP"): data[col] * data['STP_Conv'] for col in cols}
    data = data.assign(**new_data)
    return data

def Frac_conv(data, pm_col = 'PM2.5 [ug/m3 ATP]',chems= ['OC', 'EC', 'NO3', 'NH4', 'SO4'], Type = 'ATP'):
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
    cols= [c for c in data.columns.to_numpy() if (any(y in c for y in chems))&(Type in c)]
    for col in cols:
        new_col = col.replace(f' PM2.5 [ug/m3 {Type}]', "/total")
        data[new_col] = data[col]/data[pm_col]
    return data

def Cont_combo(data, cont_data):
    """
    Takes in two AQS dataframe. Combines mass concentrations
    ----------
    Paramaters
    ++++++++++
    data : [DataFrame] AQS file reorered to be indexed by datetime and possition
    cont_data : [DataFrame] PM2.5 cont dataframe

    Returns
    ++++++++++
    data : [DataFrame] AQS file with mass fraction conversions
    """
    data = data.reset_index()
    cont_data["Date(UTC)"] = pd.to_datetime(cont_data["Date(UTC)"])
    merged = data.merge(cont_data,
                        on=["Date(UTC)", "Location", 'Position'],
                        how="left")   # keeps all DF2 rows, adds PM2.5 where available
    merged = merged.set_index(["Date(UTC)", "Location", "Position"])
    merged['PM2.5 filter[ug/m3 ATP]'] = merged['PM2.5 [ug/m3 ATP]'].to_numpy()
    try:
        merged['PM2.5 [ug/m3 ATP]'] = merged['PM2.5 [ug/m3 ATP]'].combine_first(merged['PM2.5 accept[ug/m3 ATP]'])
    except:
        print('no AQS mass')
    try:
        merged['PM2.5 [ug/m3 ATP]'] = merged['PM2.5 [ug/m3 ATP]'].combine_first(merged['PM2.5 cont.[ug/m3 ATP]'])
    except:
        print('no Continuous mass')
    return merged

if __name__ == "__main__":
    # file = r"C:\Users\bensy\Documents\Research\AQS_test\daily_VOCS_2025.csv"
    # out = r'C:\Users\bensy\Documents\Research\AQS_test'
    # df, file = large_adjust(file,out)

    chem_folder = r"C:\Users\bensy\Documents\Research\AQS_Chemistry"
    pm25_folder = r"C:\Users\bensy\Documents\Research\AQS_mass"
    out_folder = r"C:\Users\bensy\Documents\Research\AQS_Processed"
    large = pd.DataFrame()
    near = pd.DataFrame()
    for root, dirs, files in os.walk(pm25_folder):
        for n in files:
            yearcsv = list(n.split('_'))[-1]
            year = yearcsv.replace('.csv', '')
            # try:
            chem_file = os.path.join(chem_folder, n)
            pm25_file = os.path.join(pm25_folder, n)
            if not os.path.exists(chem_file):
                print(f'No chemistry data for {year}')
                chem_file = ''
            if not os.path.exists(pm25_file):
                print(f'No mass data for {year}')
                pm25_file =''
            if pm25_file == chem_file:
                continue
            d,s,combo,avery = datetime_adjust(chem_file, pm25_file,out_folder)
            if large.empty:
                large = combo
            else:
                large = pd.concat([large, combo])
            if near.empty:
                near = avery
            else:
                near = pd.concat([near, avery])
            # except Exception as e: 
            #     print(f'Failed to processes files for {year} due to {e}')
    large.to_csv(join(out_folder,f"AQS_Monthly_Data_Total.csv"))
    near.to_csv(join(out_folder,f"AQS_Monthly_Data_Avery.csv"))
    print(f'Files generated at {out_folder}')