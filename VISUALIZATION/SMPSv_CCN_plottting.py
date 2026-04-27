"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
from SMPSvCCNplot_gen import line_call, hist_call,scat_call, box_call,chem_line_call,chem_scat_call, cor_scat_call, cor_box_call, cor_line_call

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
    master = master.resample(freq).apply(np.nanmean)
    master = master.dropna()
    return master,master.columns.to_numpy()

def smps_data(files,freq='d',ss = [0.1,0.7]):
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
    smps = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        file=file.set_index("DateTime Sample Start") #Set index
        if i == 0:
            smps= file
        else:
            smps = pd.concat([smps,file])
    smps.index = pd.to_datetime(smps.index)
    numsmps = [s for s in smps.columns.to_numpy() if ('.' in s) and (s.split('.')[0].isdigit())]

    # IMPORTANT: sort numerically
    numsmps = sorted(numsmps, key=lambda x: float(x))
    smps = smps[numsmps]

    # Convert to diameters 
    dp = np.array([float(n) for n in numsmps])
    logdp = np.log10(dp)
    dlogdp = np.diff(logdp)
    dlogdp = np.append(dlogdp, dlogdp[-1])  # pad last bin
    weighted = smps * dlogdp
    cols = []
    # pad to match length (assume last bin same width as previous)
    dlogdp = np.append(dlogdp, dlogdp[-1])
    for n in numsmps:
        col = f'>{float(n)}nm'
        cols.append(col)
        # bins greater than threshold
        mask = dp > float(n)
        # apply weighted sum instead of raw sum
        smps[col] = weighted.loc[:, np.array(numsmps)[mask]].sum(axis=1)
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
    cols = ['T(C)_inlet','T1(C)','T(C)_sample','T(C)_OPC','T(C)_nafion','Q(lpm)_sample','Q(lpm)_sheath','P(hPA)_sample']
    for s in ss:
        cols.append(f'N(cm-3)_avg_setpt{s}')
        cols.append(f'ss(%)_calc_setpt{s}')
        for c in ccn.columns.to_numpy(): 
            if f'N(cm-3)_cor_setpt{s}' in c:
                cols.append(c)
    ccn = ccn[cols]
    ccn = ccn.resample(freq).apply(np.nanmean)
    ccn.index.names = ['Date']
    return ccn,cols

def comb_files(smps_files,ccn_files, ss = [0.1,0.7], freq = 'd', chem = 0, cortype = ''):
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
    ccn, ccn_cols = ccn_data(ccn_files,freq,ss, cortype=cortype)
    smps, smps_cols = smps_data(smps_files,freq,ss)
    data = pd.merge(ccn[ccn_cols],smps[smps_cols],left_index = True, right_index = True)
    print(chem)
    if isinstance(chem,str):
        acsm, spec = master_data(chem)
        data = pd.merge(data, acsm ,left_index = True, right_index = True)
    return data

def plot_gen(data, mode = 0,vars = ['ss'], date = 0, group ='all', thresh = 0, cormode = False):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    vars : [list of str] columns to use while plotting (default = ['ss'])
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
    ss2nm = {'0.1':200, '0.6':100, '0.7':80}
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
    # if cormode:
    #     ccn_cors = [list(col.split('setpt'))[-1] for col in cols if "N(cm-3)_cor_" in col]
    #     ccn_cols = []
    #     choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        if choice == 'all':
            for ss in ss_vals:
                ccn_col = f'N(cm-3)_cor_setpt{ss}'
                smps_col = f'>{ss2nm[ss]}nm'
                ccn_cols.append(ccn_col)
                slct[ccn_col] = smps_col
        else:
            ccn_col = f'N(cm-3)_cor_setpt{choice}'
            smps_col = f'>{ss2nm[choice]}nm'
            slct[ccn_col] = smps_col
            ccn_cols.append(ccn_col)
    for c in ccn_cols:
        drop_indices = data[data[c]<thresh].index   
        data = data.drop(drop_indices)
    if ('Q' in vars) & (mode == 'line'): slct['Q(lpm)_sample'] = 0
    if ('T' in vars) & (mode == 'line'): slct['T(C)_sample'] = 0 
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            if cormode:
                cor_line_call(data,slct)
            else:
                line_call(data,slct)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            if cormode:
                cor_scat_call(data,slct)
            else:
                scat_call(data, slct)
        elif mode == 'hist': #histogram
            hist_call(data,slct)
        elif mode == 'box':# box and whisker
            if cormode:
                cor_box_call(data,slct)
            else: 
                box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                if cormode:
                    cor_line_call(Ydata,slct,append)
                else:
                    line_call(Ydata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                if cormode:
                    m,b,cor = cor_scat_call(Ydata, slct,append)
                else:
                    m,b,cor = scat_call(Ydata, slct,append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                hist_call(Ydata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                if cormode: 
                    cor_box_call(Ydata,slct,append)
                else:
                    box_call(Ydata,slct,append)
    elif group == 'month':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    if cormode:
                        cor_line_call(Mdata,slct,append)
                    else:
                        line_call(Mdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    if cormode:
                        m,b,cor = cor_scat_call(Mdata, slct, append)
                    else:
                        m,b,cor = scat_call(Mdata, slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    hist_call(Mdata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    if cormode:
                        cor_box_call(Mdata,slct,append)
                    else:
                        box_call(Mdata,slct,append)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    if cormode:
                        cor_line_call(Sdata,slct,append)
                    else:
                        line_call(Sdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    if cormode:
                        m,b,cor = scat_call(Sdata,slct, append)
                    else:
                        m,b,cor = scat_call(Sdata,slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    hist_call(Sdata,slct,append)
        elif mode == 'box': #box and whisker plots
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    if cormode:
                        cor_box_call(Sdata,slct,append)
                    else: 
                        box_call(Sdata,slct,append)

def chem_plot_gen(data, mode = 'line',vars = ['ss'], date = 0, group ='all', thresh = 0, chem = 'org'):
    '''
    Takes in a dataframe of SMPS and CCN and chemistry data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = line)
    vars : [list of str] columns to use while plotting (default = ['ss'])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)
    chem : [str] Species to compare data to (default = 'org')

    Returns
    ++++++++++
    none 
    '''
    ss2nm = {'0.1':200, '0.6':100, '0.7':80}
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
    comp_cols =[]
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
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        if choice == 'all':
            for ss in ss_vals:
                ccn_col = f'N(cm-3)_cor_setpt{ss}'
                smps_col = f'>{ss2nm[ss]}nm'
                new_col = f'CCN and SMPS %dev at ss={ss}'
                A = data[ccn_col].to_numpy()
                B = data[smps_col].to_numpy()
                perc_dev = (np.abs(A-B))/((A+B)/2)*100
                data[new_col] = perc_dev
                comp_cols.append(new_col)
        else:
            ccn_col = f'N(cm-3)_cor_setpt{choice}'
            smps_col = f'>{ss2nm[choice]}nm'
            new_col = f'CCN and SMPS %dev at ss={choice}'
            A = data[ccn_col].to_numpy()
            B = data[smps_col].to_numpy()
            perc_dev = (np.abs(A-B))/((A+B)/2)*100
            data[new_col] = perc_dev
            comp_cols.append(new_col)
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            chem_line_call(data,comp_cols, chem=chem)
        elif mode =='scat':
            chem_scat_call(data, comp_cols, chem=chem)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                chem_line_call(Ydata,comp_cols,append=append, chem=chem)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                m,b,cor = chem_scat_call(Ydata,comp_cols,append=append, chem=chem)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    chem_line_call(Sdata,comp_cols,append=append, chem=chem)
        elif mode == 'scat': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    chem_scat_call(Sdata,comp_cols,append=append, chem=chem)

def plot_gen_corr(data, mode = 0,vars = ['ss'], date = 0, group ='all', thresh = 0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    mode : [str] Plotting style (line,scat,hist) (default = 0, takes user input)
    vars : [list of str] columns to use while plotting (default = ['ss'])
    date : [list of str] date range for plotting (default = 0, takes user input)
            + if date = "date" or ['date'], assumed to be start date ['date':]
            + if date = ['date0','date1'], use dates contained within daterange
    group : [str] Time period to generate plots for (default = 'all')
            + 'all' - plot over whole time period given
            + 'year' - generate 1 plot per year if multiple years in data
            + 'season' - generate 1 plot per season
            + 'month' - generate 1 plot per month
    drop0s : [bool] Drop zeros in CCN data to clean(default = True)

    Returns
    ++++++++++
    none 
    '''
    ss2nm = {'0.1':200, '0.6':100, '0.7':80}
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
    smps_cols =[]
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
        mode = input('What style plot would you like to generate? (line, scat, hist, box) ')
    if 'ss' in vars: #if ss values used for plotting
        ss_vals = [list(col.split('setpt'))[-1] for col in cols if "ss(%)" in col]
        choice = input(f'Which ss% value would you like to use? ({', '.join(ss_vals)}, or all) ')
        if choice == 'all':
            for ss in ss_vals:
                ccn_cols = [col for col in data.columns.to_numpy() if f'N(cm-3)_cor_setpt{ss}' in col]
                smps_col = f'>{ss2nm[ss]}nm'
                smps_cols.append(smps_col)
                slct[smps_col] = ccn_cols
        else:
            ccn_cols = [col for col in data.columns.to_numpy() if f'N(cm-3)_cor_setpt{choice}' in col]
            smps_col = f'>{ss2nm[choice]}nm'
            smps_cols.append(smps_col)
            slct[smps_col] = ccn_cols
    for c in ccn_cols:
        drop_indices = data[data[c]<thresh].index   
        data = data.drop(drop_indices)
    if ('Q' in vars) & (mode == 'line'): slct['Q(lpm)_sample'] = 0
    if ('T' in vars) & (mode == 'line'): slct['T(C)_sample'] = 0 
    if group =='all':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            line_call(data,slct)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            scat_call(data, slct)
        elif mode == 'hist': #histogram
            hist_call(data,slct)
        elif mode == 'box':# box and whisker
            box_call(data,slct)
    elif group == 'year':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                line_call(Ydata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                m,b,cor = scat_call(Ydata, slct,append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                hist_call(Ydata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                append=f"{year}"
                box_call(Ydata,slct,append)
    elif group == 'month':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    line_call(Mdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    m,b,cor = scat_call(Mdata, slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    hist_call(Mdata,slct,append)
        elif mode == 'box': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for month in Ydata.month.unique():
                    Mdata=Ydata[Ydata.month == month]
                    append=f"{month}/{year}"
                    box_call(Mdata,slct,append)
    elif group == 'season':
        if mode == 'line': ## Line plot SMPS and CCN vs date
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    line_call(Sdata,slct,append)
        elif mode == 'scat': # scatter plot SMPS vs CCN
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season == season]
                    append=f"{season} {year}"
                    m,b,cor = scat_call(Sdata,slct, append)
        elif mode == 'hist': #histogram
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    hist_call(Sdata,slct,append)
        elif mode == 'box': #box and whisker plots
            for year in data.year.unique():
                Ydata=data[data.year == year]
                for season in Ydata.season.unique():
                    Sdata=Ydata[Ydata.season== season]
                    append=f"{season} {year}"
                    box_call(Sdata,slct,append)

if __name__ == '__main__':
    smps = [r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv"]#list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv"]#list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    master =  r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"
    data = comb_files(smps,ccn, freq='d')
    data_chem = comb_files(smps,ccn, chem= master, freq="d")
    out = r"C:\Users\bensy\Documents\Research\CCN_vSMPS.csv"#input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        data.to_csv(out)
    # bad_dates = ['8/15/2025 00:00:00','12/01/2025 00:00:00']
    plot_gen(data,mode = 'line',thresh=5,group='year')
    # chem_plot_gen(data_chem,mode = 'scat',chem = 'org')