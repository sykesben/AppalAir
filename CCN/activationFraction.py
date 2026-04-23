"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr 
from scipy import stats as st 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as LSfit
pd.set_option('mode.chained_assignment', None)
from plotgen import box_call, line_call, hist_call,scat_call
large_nm = 60

def master_data(f,freq='d'):
    '''
    Takes in the master file and specifically cuts out the AQS data
    ----------

    Parameters
    ++++++++++
    f : [list of str] Paths to Master file
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    master : [DataFrame] Master data file
    spec : [list of str] Names of used columns from chemistry output
    '''
    master=pd.read_csv(f) #read in AQS file
    master=master.set_index("Local time (UTC-5)") #Set index
    master['Date(UTC)'] = pd.to_datetime(master.index) + pd.Timedelta(hours=5)
    specs = ['NH4_11000','SO4_11000','NO3_11000','Org_11000','1hrMC_µg/m3','org/total','SO4/total']
    master=master.set_index("Date(UTC)")
    master = master[specs]
    master.columns = master.columns.str.replace('_11000', ' [µg/m3] ACSM')
    master = master.resample(freq).mean()
    master = master.dropna()
    return master,master.columns.to_numpy()

def smps_data(files,freq='d',nm = [large_nm]):
    '''
    Takes in a list of smps files and filters out the particle size concentration depending
    on comparable ss% values.
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame
    ss : [list of floats] ss% set points from CCN

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    global large_nm
    smps = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        file=file.set_index("DateTime Sample Start") #Set index
        if i == 0:
            smps= file
        else:
            smps = pd.concat([smps,file])
    numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]

    # IMPORTANT: sort numerically
    numsmps = sorted(numsmps, key=lambda x: float(x))
    smps_nums = smps[numsmps]

    # Convert to diameters 
    dp = np.array([float(n) for n in numsmps])
    logdp = np.log10(dp)
    dlogdp = np.diff(logdp)
    dlogdp = np.append(dlogdp, dlogdp[-1])  # pad last bin
    weighted = smps_nums * dlogdp
    cols = ['Median (nm)','Mean (nm)','Geo. Mean (nm)','Mode (nm)','Geo. Std. Dev','Total Concentration (#/cm³)']
    # pad to match length (assume last bin same width as previous)
    dlogdp = np.append(dlogdp, dlogdp[-1])
    for n in nm:
        col = f'>{float(n)}nm'
        percol = f'>{float(n)}nm[%]'
        cols.append(col)
        cols.append(percol)
        # bins greater than threshold
        mask = dp > float(n)
        # apply weighted sum instead of raw sum
        wSum = weighted.loc[:, np.array(numsmps)[mask]].sum(axis=1)
        smps[col] = wSum
        smps[percol] = wSum/smps['Total Concentration (#/cm³)']*100

    smps.index = pd.to_datetime(smps.index)
    smps = smps[cols]
    smps = smps.resample(freq).mean()
    smps.index.names = ['Date']
    return smps,cols

def ccn_data(files, freq ='d', ss = [0.1,0.7], cortype =''):
    '''
    Takes in a list of CCN files and returns a processed dataframe with important columns 
    for plotting or further analysis
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrame (default = 'd')
    ss : [list of floats] ss% set points from CCN (default = [0.1,0.7])
    cortype : [float] if specific correction column is wanted (default = '')

    Returns
    ++++++++++
    ccn : [DataFrame] Combined CCN data from all inputted files
    cols : [list of str] Names of used columns from CCN output
    '''
    ccn = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in ccn file
        try:
            file=file.set_index('Datetime(UTC)') #Set index
        except:
            try:
                file=file.set_index('Datetime UTC') #Set index
            except:
                file = file.set_index('Date String (YYYY-MM-DD hh:mm:ss) UTC')
        file.index = file.index.rename('Datetime(UTC)')
        if i == 0:
            ccn = file
        else:
            ccn = pd.concat([ccn,file])
    ccn.index = pd.to_datetime(ccn.index, format='mixed')
    cols = []
    ss_cols = []
    for c in ccn.columns.to_numpy(): 
        if (f'N(cm-3)_cor_setpt' in c) | (f'ss(%)_calc_setpt' in c) | (f'N(cm-3)_avg_setpt' in c):
            cols.append(c)
        if (f'N(cm-3)_cor_setpt' in c):
            ss_cols.append(c)
    ccn = ccn[cols]
    ccn = ccn.resample(freq).mean()
    ccn.index.names = ['Date']
    return ccn,cols,ss_cols

def comb_files(smps_files,ccn_files, freq = 'd', chem = 0, cortype = ''):
    '''
    Takes in a list of CCN and SMPS files and returns a combined dataframe with important columns 
    from both for plotting or further analysis
    ----------

    Parameters
    ++++++++++
    smps_files : [list of str] Paths to SMPS files
    ccn_files : [list of str] Paths to CCN files
    freq : [str] Resample frequency for DataFrames (default = 'd')
    ss : [list of floats] ss% set points from CCN (default = [0.1,0.7])
    chem : [list of str] Paths to Chemistry data if included (default = 0)
    cortype : [float] if specific correction column is wanted (default = '')

    Returns
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data from all inputted files
    '''
    ccn, ccn_cols,ss_cols = ccn_data(ccn_files,freq, cortype=cortype)
    smps, smps_cols = smps_data(smps_files,freq)
    data = pd.merge(ccn[ccn_cols],smps[smps_cols],left_index = True, right_index = True)
    print(chem)
    if isinstance(chem,str):
        acsm, spec = master_data(chem)
        data = pd.merge(data, acsm ,left_index = True, right_index = True)
    return data,ss_cols,smps_cols

def find_activation(data, smps_cols, ss_cols, ss_vals = ['0.1','0.15','0.4','0.7']):
    FA_cols= []
    for col in ss_cols:
        print(col)
        ss = col.split('cor_setpt')[-1] 
        if ss in ss_vals:
            cnt_at_ss = data[col].to_numpy()
            print(col)
            # print(cnt_at_ss)
            cnt_tot= data['Total Concentration (#/cm³)'].to_numpy()
            mask = cnt_tot != 0 
            Facts = np.divide(cnt_at_ss, cnt_tot, where=cnt_tot!=0, out=np.zeros_like(cnt_tot, dtype=float))
            Facts = np.array(Facts)
            data[f'Fact_at_{ss}'] = Facts
            FA_cols.append(f'Fact_at_{ss}')
    return data, FA_cols

def cov(x):
    x = x.dropna()
    mean = x.mean()
    if mean == 0:
        return np.nan  # avoid divide-by-zero
    return x.std() / mean*100

def find_deviation(data, ss_cols,const_room =0.10, thresh =0.1, const = 'GeoMean',Var= 'Org', ss_vals = ['0.1','0.15','0.7']):
    FA_cols= []
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        if ss in ss_vals:
            Fact = f'Fact_at_{ss}'
            FA_cols.append(Fact)
            title = ['Y',' vs ','X']
            y_label = 'F_act CoV[%]'
            if const =='Org':
                col = 'org/total'
                title[0] = 'Fact CoV at avg Forg'
            elif const =='SO4':
                col = 'SO4/total'
                title[0] = 'Fact CoV at avg Fso4'
            elif const == 'GeoMean':
                col = f'Geo. Mean (nm)'
                title[0] = 'Fact CoV at avg Geo. Mean'
            if Var == 'Org':
                var = 'org/total'
                x_label = 'F_org[1]'
                title[-1] = 'Organic Fraction'
            elif Var == 'SO4':
                var = 'SO4/total'
                x_label = 'F_so4[1]'
                title[-1] = 'Sulfate Fraction'
            elif Var == 'GeoMean':
                var = f'Geo. Mean (nm)'
                x_label= 'Geo Mean [nm]'
                title[-1] = 'Geometric Mean'
            med = data[col].median()
            mask = (data[col] - med).abs() <= (abs(med) * const_room)
            data = data[mask]
            print(data[var].min())
            data['decile'] = pd.qcut(data[var], 10, labels=False)
            result = data.groupby('decile').agg({
                        Fact: cov,
                        var: 'mean'})
            drop_indices = data[data[Fact]<thresh].index   
            data = data.drop(drop_indices)
            print(result)
            slct ={}
            slct[Fact] = var
            scat_call(result, slct, x_label=x_label ,y_label=y_label,title=''.join(title), verbose=True,single=True)
    return data, FA_cols


def plot_gen(data, mode = 0,yvars = 'Fa',xvars = 'Dt',const = '',const_room =0.05,ss_vals = ['0.7'], lines ='',date = 0, group ='all', thresh = 0, cormode = False,showCor=False, singleLine=True):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    yvars : [str] y vars to use while plotting (default = 'Fa')
    xvars : [str] x vars to use while plotting (default = 'Dt'/Datetime)
    const : [str] Get values only where a constant is held (defualt =''/none)
    ss : [list of float] list of ss values to use (default = [0.1,0.4,0.7])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)
    cormode : [bool] for comparing between corrections (default = False)

    Returns
    ++++++++++
    none 
    '''
    global large_nm
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    data['year'] = data.index.year.to_numpy()
    data['month']= data.index.month.to_numpy()
    data['season'] = [seasons(n) for n in data.index.month.to_numpy()]
    slct = {}
    cols = data.columns.to_numpy()
    ccn_cols =[]
    if date != 0: #if date is passed, split the data using the passed date range
        if isinstance(date, list):
            if len(date) > 1:
                date0 = pd.to_datetime(date[0])
                date1 = pd.to_datetime(date[1])
                data = data.loc[date0:date1]
            else:
                date0 = pd.to_datetime(date[0])
                data = data.loc[date0:]
        elif isinstance(date, str):
            date0 = pd.to_datetime(date)
            data = data.loc[date0:]
    if mode == 0: #if default value used and no mode passed, user input mode
        mode = input('What style plot would you like to generate? (line, scat, hist)')
    y_cols = []
    x_cols = []
    split = 0 
    y_name = ''
    x_name = ''
    title_list = ['Y', ' vs ', 'X']
    '''Lines to split on'''
    if lines == 'Org':
        sorted_array = np.sort(data['org/total'].to_numpy())
        N25 = int(len(sorted_array) * 0.25) #lower 25% of org%
        N75 = int(len(sorted_array) * 0.75) #top 25% of org%
        O25 = sorted_array[:N25]
        O50 = sorted_array[N25:N75]
        O75 = sorted_array[N75:]
        split = {}
        split['org/total'] = ['Bottom 25% of Org',O25], ['Top 25% of Org', O75]
    elif lines == 'GeoMean':
        sorted_array = np.sort(data[f'Geo. Mean (nm)'].to_numpy())
        N25 = int(len(sorted_array) * 0.25) #lower 25% of GeoMean
        N75 = int(len(sorted_array) * 0.75) #top 25% of GeoMean
        O25 = sorted_array[:N25]
        O50 = sorted_array[N25:N75]
        O75 = sorted_array[N75:]
        split = {}
        split[f'Geo. Mean (nm)'] = ['Bottom 25% of GeoMean',O25], ['Top 25% of GeoMean', O75]
    # input(list(split.values()))
    '''Set Y values'''
    if yvars =='Fa': #if activation fraction values used for plotting
        y_name = 'Fa [1]'
        title_list[0] = ('Activation Fraction')
        if len(ss_vals) > 1:
            choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        else: choice = ss_vals[0]
        if choice == 'all':
            for ss in ss_vals:
                Fa = f'Fact_at_{ss}'
                y_cols.append(Fa)
        else:
            Fa = f'Fact_at_{choice}'
            y_cols.append(Fa)
    elif yvars =='Std_Fa': #if activation fraction standard deviation used for plotting
        y_name = 'Fa Dev[1]'
        title_list[0] = ('Fact St.Dev.')
        if len(ss_vals) > 1:
            choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        else: choice = ss_vals[0]
        if choice == 'all':
            for ss in ss_vals:
                Fa = f'Fact_at_{ss}'
                y_cols.append(Fa)
        else:
            Fa = f'Fact_at_{choice}'
            y_cols.append(Fa)
    for c in y_cols:
        drop_indices = data[data[c]<thresh].index   
        data = data.drop(drop_indices)
    '''Set constants'''
    if const =='Org':
        #only keep values where org fraction is close to median value
        col ='org/total'
        med = data[col].median()
        mask = (data[col] - med).abs() <= (abs(med) * const_room)
        data = data[mask]
        y_name += ' at const org%'
    elif const == 'GeoMean':
        #only keep values where org fraction is close to median value
        col = f'Geo. Mean (nm)'
        med = data[col].median()
        mask = (data[col] - med).abs() <= (abs(med) * const_room)
        data = data[mask]
        y_name += ' at const Geo. Mean'

    '''Set X values'''
    if xvars=='Fa': #if activation fraction values used for plotting
        x_name = 'Fa [1]'
        title_list[-1] = ('Activation Fraction')
        if len(ss_vals) > 1:
            choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        else: choice = ss_vals[0]
        if choice == 'all':
            for ss in ss_vals:
                Fa = f'Fact_at_{ss}'
                x_cols.append(Fa)
        else:
            Fa = f'Fact_at_{choice}'
            x_cols.append(Fa)
    elif xvars=='Org': #if organic fraction values used for plotting
        x_name = 'Org [1]'
        title_list[-1] = ('Organic Fraction')
        org = f'org/total'
        for i in range(len(y_cols)):
            x_cols.append(org)
    elif xvars=='SO4': #if organic fraction values used for plotting
        x_name = 'SO4 [1]'
        title_list[-1] = ('Sulfate Fraction')
        sulf = f'SO4/total'
        for i in range(len(y_cols)):
            x_cols.append(sulf)
    elif xvars=='NO3': #if Nitrate fraction values used for plotting
        x_name = 'NO3 [1]'
        title_list[-1] = ('Nitrate Fraction')
        sulf = f'NO3/total'
        for i in range(len(y_cols)):
            x_cols.append(sulf)
    elif xvars=='Large': #if percent above threshold used for plotting
        x_name = f'>{large_nm}nm [%]'
        title_list[-1] = ('Percent of Large Particles')
        size = f'>{float(large_nm)}nm[%]'
        for i in range(len(y_cols)):
            x_cols.append(size)
    elif xvars=='GeoMean': #if geometric mean values used for plotting
        x_name = 'Geo. Mean [nm]'
        title_list[-1] = ('Geometric Mean')
        GM = f'Geo. Mean (nm)'
        for i in range(len(y_cols)):
            x_cols.append(GM)
    elif xvars == 'Dt': #if Datetime value used for plotting
        title_list[1] =''
        title_list[-1] =''
        x_name = 'Datetime'
        x_cols = (np.zeros(len(y_cols))).astype(int)
    title = f'{''.join(title_list)}'
    if (mode == 'box')|(mode=='scat'):
        for i in range(len(y_cols)):
            slct[y_cols[i]] = x_cols[i]
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            line_call(data,slct,y_label=y_name, x_label=x_name, title= title)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            scat_call(data, slct, y_label=y_name, x_label=x_name, title= title,split=split,verbose=showCor, single=singleLine)
        elif mode == 'box':# box and whisker
            box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f" {year}"
                title_new = title + append
                line_call(Ydata,slct,y_label=y_name, x_label=x_name, title= title_new)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f" {year}"
                title_new = title + append
                m,b,cor = scat_call(Ydata, slct, y_label=y_name, x_label=x_name, title= title_new, split=split,verbose=showCor, single=singleLine)
                Ydata=data[data.year == year]
                append=f"{year}"
                hist_call(Ydata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f" {year}"
                title_new = title + append
                box_call(Ydata,slct,y_label=y_name,title= title_new)
    elif group == 'month':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f" {month}/{year}"
                    title_new= title +append
                    line_call(Mdata,slct,y_label=y_name, x_label=x_name, title= title_new)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f" {month}/{year}"
                    title_new = title + append
                    m,b,cor = scat_call(Mdata, slct, y_label=y_name, x_label=x_name, title= title_new, split=split,verbose=showCor, single=singleLine)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f" {month}/{year}"
                    title_new = title + append
                    box_call(Mdata,slct,y_label=y_name,title= title_new)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f" {season} {year}"
                    title_new = title + append
                    line_call(Sdata,slct, y_label=y_name, x_label=x_name, title= title_new)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    # input(Sdata.columns)
                    append=f" {season} {year}"
                    title_new = title + append
                    m,b,cor = scat_call(Sdata,slct, y_label=y_name, x_label=x_name, title= title_new, split=split,verbose=showCor, single=singleLine)
        elif mode == 'box': #box and whisker plots
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f" {season} {year}"
                    title_new = title + append
                    box_call(Sdata,slct, y_label=y_name,title= title_new)

if __name__ == '__main__':
    bad_dates = [[pd.to_datetime('10/01/2025 00:00:00'),pd.to_datetime('12/15/2025 00:00:00')]]
    smps = [r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv"]#list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv"]#list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    master =  r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"
    data,ss_cols,smps_cols = comb_files(smps,ccn, freq='d')#chem=master
    mask = pd.Series(False, index=data.index)
    for date in bad_dates:
        mask |= (data.index >= date[0]) & (data.index <= date[-1])
    data = data[~mask]
    data, FA_cols = find_activation(data, smps_cols, ss_cols)
    # data, FA_cols = find_deviation(data, ss_cols)
    plot_gen(data,mode='scat',xvars='GeoMean', group='all',thresh=0.10, ss_vals=['0.7'],showCor=True,singleLine =True) # 
    input(data)
    out = r"C:\Users\bensy\Documents\Research\CCN_activation_Fraction.csv"#input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        data.to_csv(out)

