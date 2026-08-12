"""
Date: 3/3/26
Author: Ben Sykes
Purpose: generate plots between CCN and SMPS
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
from scipy.stats import linregress, pearsonr 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as LSfit
pd.set_option('mode.chained_assignment', None)
plt.rcParams['font.size'] = 12

def critical_diameter(ss, kappa=0.1, T=298):
    """
    Estimated critical diameter in nm from SS (%)
    """
    sigma = 0.072  # surface tension (N/m)
    Mw = 0.018     # kg/mol
    R = 8.314
    rho_w = 1000   # kg/m3

    A = (4 * sigma * Mw) / (R * T * rho_w)
    ss = float(ss) / 100  # % to fraction
    Dcrit = ((4 * A**3) / (27 * kappa * (np.log(1 + ss))**2))**(1/3)
    return Dcrit * 1e9  # m to nm

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

def cpc_data(f, freq = 'd'):
    file =pd.read_csv(f, skiprows=1)
    file=file.set_index('DateTimeUTC')
    file.index = pd.to_datetime(file.index, format='mixed')
    file = file.resample(freq).mean()
    file.index.names = ['Date']
    file.rename(columns={'N_N71': 'CPC N[cm-3]'}, inplace=True)
    return file


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
    total = smps['Total Concentration (#/cm³)'].to_numpy()
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
    for n in numsmps:
        col = f'>{float(n)}nm'
        cols.append(col)
        # bins greater than threshold
        mask = dp >= float(n)
        # apply weighted sum instead of raw sum
        smps[col] = weighted.loc[:, np.array(numsmps)[mask]].sum(axis=1)
    smps['Total Concentration (#/cm³)'] = total
    cols.append('Total Concentration (#/cm³)')
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

def comb_files(smps_files,ccn_files, cpc_file, freq = 'd', chem = 0, cortype = '',extra_fout = 0,ss_vals = ['0.1','0.15','0.25','0.4','0.7']):
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
    smps, diam_cols = smps_data(smps_files,freq)
    cpc = cpc_data(cpc_file, freq)
    input(smps)
    data = pd.merge(ccn[ccn_cols],smps[diam_cols],left_index = True, right_index = True)
    data = pd.merge(data,cpc,left_index = True, right_index = True)
    for ss in ss_vals:
        data[f'F_act_SMPS(ss={ss}%)'] = data[f'N(cm-3)_cor_setpt{ss}']/data[f'Total Concentration (#/cm³)']
        data[f'F_act({ss}%)'] = data[f'N(cm-3)_cor_setpt{ss}']/data[f'CPC N[cm-3]']
    print(chem)
    if isinstance(chem,str):
        acsm, spec = master_data(chem)
        data = pd.merge(data, acsm ,left_index = True, right_index = True)
    if extra_fout != 0:
        data.to_csv(extra_fout)
    return data,ss_cols,diam_cols


def find_act(data, diams, ss_cols, key ='', Kappa = 0.1):
    cut_off = {}
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        est_D = critical_diameter(ss,kappa=Kappa) # estimated critical Diameter
        Fact_at_ss = data[f'F_act{key}({ss}%)']
        devs = []
        for diam in diams:
            cnt_at_diam = data[diam]
            dev = cnt_at_ss-cnt_at_diam*.50 #50% of the total particles can activate as CCN
            devs.append(float(dev))
        devs = np.array(devs)
        devs = devs[:-1]
        indices = np.where(~np.isnan(devs))[0]  # Get integer indices
        devs = devs[indices]  # Use integer indices to filter 

        notFound = True
        start = 0
        end = len(devs)-1
        while notFound: 
            mid = (start + end)//2
            s_dev = devs[start]
            e_dev = devs[end]
            m_dev = devs[mid]
            if s_dev >0:
                notFound = False
                diam_top = float(diams[start].replace('>','').replace('nm',''))
                diam_bot = 0#float(diams[start].replace('>','').replace('nm',''))
                diam_mid = find_mid(diam_top,diam_bot, float(devs[start]),float(devs[start]))
                cut_off[ss] = [diam_mid, diam_bot,diam_top]
            # if ss == '0.7':
                # print(ss)
                # print(f'diameters = {diams[start]}-{diams[end]}')
                # input(f'dev = {devs[start]} , {devs[mid]} , {devs[end]}')
            if  m_dev * s_dev > 0: #same sign
                start = mid
            elif m_dev * e_dev > 0: #same sign
                end = mid
            if abs(start-end) ==1:
                notFound = False
                diam_top = float(diams[end].replace('>','').replace('nm',''))
                diam_bot = float(diams[start].replace('>','').replace('nm',''))
                diam_mid = find_mid(diam_top,diam_bot, float(devs[end]),float(devs[start]))
                cut_off[ss] = [diam_mid, diam_bot,diam_top]
    return cut_off

def find_mid(D_top, D_bot, diff_top, diff_bot):
    delta_D = abs(D_top-D_bot)
    delta_diff = abs(diff_bot) + abs(diff_top)
    dev_from_bottom = abs(diff_bot)/delta_diff
    diameter = delta_D*dev_from_bottom + D_bot
    return diameter

def find_activation(data, diams, ss_cols):
    act_perc = {}
    d_index = [float(x.replace('>','').replace('nm','')) for x in diams[:-1]]
    for col in ss_cols:
        ss = col.split('cor_setpt')[-1] 
        # est_D = critical_diameter(ss,kappa=Kappa) # estimated critical Diameter
        Acts = data[f'F_act({ss}%)'].to_numpy()
        Acts = np.array(Acts)
        Acts = Acts[:-1]
        indices = np.where(~np.isnan(Acts))[0]  # Get integer indices
        Acts = Acts[indices]  # Use integer indices to filter 
        act_perc[f'{ss}'] = Acts
    return act_perc

def cut_off_curve(data,freq = 'M', group ='all'):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
    group : [str] Time period to generate plots for (default = 'all')
        + 'all' - plot over whole time period given
        + 'year' - generate 1 plot per year if multiple years in data
        + 'season' - generate 1 plot per season
        + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    df = pd.DataFrame()
    if freq =='S':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        month_vals = np.array([(c.split('_')[-1]).split('/')[0] for c in cols])
        szns= [seasons(float(c)) for c in month_vals]
        Scols = [' - '.join([seasons(float(x)), y]) for x, y in zip(month_vals,year_vals)]
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        szn_cols = []
        for yr in np.unique(year_vals):
            for szn in np.unique(szns):
                s_idx = [i for i,c in enumerate(Scols) if (szn in c)&(yr in c)]
                # print(Scols)
                # input(s_idx)
                cmn_cut = list(set(s_idx)&set(cutoffs))
                cmn_low = list(set(s_idx)&set(lowers))
                cmn_up = list(set(s_idx)&set(uppers))
                cut_col = [cols[i] for i in cmn_cut]
                low_col = [cols[i] for i in cmn_low]
                up_col = [cols[i] for i in cmn_up]    
                data[f'D_cutoff_{szn}/{yr}'] = data[cut_col].mean(axis =1)
                data[f'D_lower_{szn}/{yr}'] = data[low_col].mean(axis =1)
                data[f'D_upper_{szn}/{yr}'] = data[up_col].mean(axis =1)
                szn_cols.extend([f'D_cutoff_{szn}/{yr}',f'D_lower_{szn}/{yr}',f'D_upper_{szn}/{yr}'])
        data = data[szn_cols]
    if freq =='Y':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        yr_cols = []
        for yr in np.unique(year_vals):
            y_idx = [i for i,c in enumerate(year_vals) if yr in c]
            cmn_cut = list(set(y_idx)&set(cutoffs))
            cmn_low = list(set(y_idx)&set(lowers))
            cmn_up = list(set(y_idx)&set(uppers))
            cut_col = [cols[i] for i in cmn_cut]
            low_col = [cols[i] for i in cmn_low]
            up_col = [cols[i] for i in cmn_up]            
            data[f'D_cutoff_{yr}'] = data[cut_col].mean(axis =1)
            data[f'D_lower_{yr}'] = data[low_col].mean(axis =1)
            data[f'D_upper_{yr}'] = data[up_col].mean(axis =1)
            yr_cols.extend([f'D_cutoff_{yr}',f'D_lower_{yr}',f'D_upper_{yr}'])
        data = data[yr_cols]
    bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
    slct = [c for c in data.columns.to_numpy() if ('cutoff' in c)&(not any(y in c for y in bad_szns))]# & ('autumn' not in c)] 
    scat_call(data, slct)
    input(data)

def Fact_curve(data,freq = 'M', group ='all', ss_choices = ['0.1','0.4','0.7']):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
    group : [str] Time period to generate plots for (default = 'all')
        + 'all' - plot over whole time period given
        + 'year' - generate 1 plot per year if multiple years in data
        + 'season' - generate 1 plot per season
        + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    if len(ss_choices)>1:
        for ss in ss_choices:
            if freq =='S':
                year_vals = np.array([c.split('/')[-1] for c in cols])
                month_vals = np.array([(c.split('_')[-1]).split('/')[0] for c in cols])
                ss_vals = np.array([(c.split('_')[-2]) for c in cols])
                szns= [seasons(float(c)) for c in month_vals]
                select = ['_'.join([s,'/'.join([seasons(float(m)), y])])for s,m, y in zip(ss_vals, month_vals,year_vals)]
                szn_cols = []
                for yr in np.unique(year_vals):
                    for szn in np.unique(szns):
                        s_idx = [i for i,c in enumerate(select) if (szn in c)&(yr in c)&(c.split('_')[0] == ss)]
                        act_col = [cols[i] for i in s_idx]
                        print(act_col)
                        data[f'Fact_{szn}/{yr}'] = data[act_col].mean(axis =1)
                        szn_cols.extend([f'Fact_{szn}/{yr}'])
                plot_data = data[szn_cols]
            if freq =='Y':
                year_vals = np.array([c.split('/')[-1] for c in cols])
                cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
                lowers = [i for i,c in enumerate(cols) if 'lower' in c]
                uppers = [i for i,c in enumerate(cols) if 'upper' in c]
                yr_cols = []
                for yr in np.unique(year_vals):
                    y_idx = [i for i,c in enumerate(year_vals) if (yr in c)&(c.split('_')[0] == ss)]
                    act_col = [cols[i] for i in y_idx]
                    data[f'Fact_{yr}'] = data[act_col].mean(axis =1)
                    yr_cols.extend([f'Fact_{yr}'])
                plot_data = data[yr_cols]
            bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
            slct = [c for c in data.columns.to_numpy() if ('Fact' in c) & (not any(y in c for y in bad_szns))]# & ('autumn' not in c)] 
            plot_data = plot_data[slct]
            plot_data[plot_data>2.1] = np.nan
            plot_data = plot_data.dropna(how='any')
            # input(plot_data)
            scat_call_Fact(plot_data, slct, name= 'Fact',append=ss)
    input(data)

def Dist_curve(data,freq = 'S', group ='all'):
    '''
    Takes in a dataframe of SMPS and CCN data and generates interactive plots based on 
    the chosen columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    freq : [str] How the data is averaged (default = M)
        + 'M' - month, no averaging
        + 'S' - Seasonal averaging
        + 'Y' - Yearly averaging
    group : [str] Time period to generate plots for (default = 'all')
        + 'all' - plot over whole time period given
        + 'year' - generate 1 plot per year if multiple years in data
        + 'season' - generate 1 plot per season
        + 'month' - generate 1 plot per month

    Returns
    ++++++++++
    none 
    '''
    def seasons(num):
        if num in [12,1,2]: return 'winter'
        elif num in [3,4,5]: return 'spring'
        elif num in [6,7,8]: return 'summer'
        else: return 'autumn' 
    cols = data.columns.to_numpy()
    if freq =='S':
        dt = pd.to_datetime(data.index.to_numpy())
        # input(dt.year)
        year = [str(yr).replace('20','') for yr in dt.year]
        month = dt.month
        data['season'] = ['/'.join([seasons(float(m)), y])for m, y in zip(month,year)]
        data = data.set_index('season')
        cols_select = [col for col in data.columns.to_numpy() if '>' in col]
        tot = data['Total Concentration (#/cm³)'].to_numpy()
        data = data[cols_select]
        grouped = data.groupby(data.index).mean()
        grouped.columns = [float(str(c).replace('>','').replace('nm','')) for c in grouped.columns.to_numpy()]
        input(grouped)
    if freq =='Y':
        year_vals = np.array([c.split('/')[-1] for c in cols])
        cutoffs = [i for i,c in enumerate(cols) if 'cutoff' in c]
        lowers = [i for i,c in enumerate(cols) if 'lower' in c]
        uppers = [i for i,c in enumerate(cols) if 'upper' in c]
        yr_cols = []
        for yr in np.unique(year_vals):
            y_idx = [i for i,c in enumerate(year_vals) if (yr in c)&(c.split('_')[0] == ss)]
            act_col = [cols[i] for i in y_idx]
            data[f'Fact_{yr}'] = data[act_col].mean(axis =1)
            yr_cols.extend([f'Fact_{yr}'])
        plot_data = data[yr_cols]
    bad_szns = ['spring/24','autumn/25','summer/26', 'autumn/26']
    slct = [c for c in data.columns.to_numpy() if ('Fact' in c) & (not any(y in c for y in bad_szns))]# & ('autumn' not in c)] 
    plot_data = plot_data[slct]
    plot_data[plot_data>2.1] = np.nan
    plot_data = plot_data.dropna(how='any')
    # input(plot_data)
    scat_call_Fact(plot_data, slct, name= 'Fact',append=ss)
    input(data)

def scat_call(data,slct,name = 'D50',append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    yerr = []
    fits = []
    year = 0
    clr_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for cutoff in list(slct):
        date = cutoff.split('cutoff_')[-1]
        szn, yr = date.split('/')
        color = colors[szn][yr]
        clr_out.append(color)
        up = cutoff.replace('cutoff','upper')
        low = cutoff.replace('cutoff','lower')
        upper = abs(data[up].to_numpy()-data[cutoff].to_numpy())
        lower = abs(data[low].to_numpy()-data[cutoff].to_numpy())
        label = cutoff.split('cuttoff_')[-1]
        asym = [lower,upper]
        y.append(data[cutoff].to_numpy())
        x.append(data.index.to_numpy().astype(float))
        coeff = np.polyfit(data.index.to_numpy().astype(float), data[cutoff].to_numpy(), 5)
        p = np.poly1d(coeff)
        fit = (p(data.index.to_numpy().astype(float)))
        fits.append(fit)
        yerr.append(asym)
        leg.append(f'{name} {label.split('cutoff_')[-1]}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
    scat_plot(x,y,yerr, fits, leg,clr_out, append=append)
    
def scat_plot(x,y, err,fits,legs,colors,append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    points = []
    for i in range(len(y)):
        ps,caps,bar = ax.errorbar(x[i],y[i], yerr= err[i],ls = '', marker = '.', capsize=3,color= colors[i])
        L = ax.plot(x[i], fits[i],label = legs[i], color = ps.get_color())
        lines.append(L)
        pnt = [ps]
        pnt.extend(list(caps))
        pnt.append(list(bar)[0])
        points.append(pnt)
    leg = ax.legend()
    lined = dict()
    for legline, origline,point in zip(leg.get_lines(), lines, points):       
        legline.set_picker(10)  # 10 pts tolerance
        lined[legline] = origline[0], point
        # input(lined[legline])

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline, point = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        for p in point:
            p.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel('Diameter [nm]')
    ax.set_xlabel('ss [%]')
    if append ==0:
        ax.set_title(f"Activation diameter of CCN")
    else:
        ax.set_title(f"Activation diameter of CCN for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

def scat_call_Fact(data,slct,name = 'Fact',append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    yerr = []
    fits = []
    year = 0
    clr_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for fact in list(slct):
        date = fact.split('_')[-1]
        szn, yr = date.split('/')
        color = colors[szn][yr]
        clr_out.append(color)
        label = fact.split('cuttoff_')[-1]
        y.append(data[fact].to_numpy())
        x.append(data.index.to_numpy().astype(float))
        coeff = np.polyfit(data.index.to_numpy().astype(float), data[fact].to_numpy(), 5)
        p = np.poly1d(coeff)
        fit = (p(data.index.to_numpy().astype(float)))
        fits.append(fit)
        leg.append(f'{name} {label.split('cutoff_')[-1]}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
    scat_plot_Fact(x,y, fits, leg,clr_out, append=append)

def scat_plot_Fact(x,y,fits,legs,colors,append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    points = []
    ax.hlines(y = [1,0.5], colors=['k', '#808080'], linestyles=['dashed'],xmax=np.max(x),xmin = np.min(x))
    for i in range(len(y)):
        S = ax.plot(x[i],y[i],ls = '', marker = '.',color= colors[i])
        L = ax.plot(x[i], fits[i],label = legs[i], color = colors[i])
        lines.append(L)
        points.append(S)
    leg = ax.legend()
    lined = dict()
    for legline, origline,point in zip(leg.get_lines(), lines, points):       
        legline.set_picker(10)  # 10 pts tolerance
        lined[legline] = origline[0], point
        # input(lined[legline])

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline, point = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        for p in point:
            p.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel('Fact [1]')
    ax.set_xlabel('Diameter [nm]')
    if append ==0:
        ax.set_title(f"Activation Faction of CCN")
    else:
        ax.set_title(f"Activation Fraction of CCN for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

def scat_call_dist(data,append=0):
    '''
    Takes in a dataframe of SMPS and CCN data and generates an interactive scatter plot based on 
    the selected columns and mode.
    ----------

    Parameters
    ++++++++++
    data : [DataFrame] Combined CCN and SMPS data
    slct : [dict of str] relative CCN and SMPS column names for processing
    append : [any] value to append to end of plot title (default = 0[no appending])

    Returns
    ++++++++++
    m_all : [list of float] slopes of lines
    b_all : [list of float] intercept of lines
    r_all : [list of float] pearson coefficient of lines
    '''
    x = []
    y = []
    leg = []
    yerr = []
    fits = []
    year = 0
    clr_out = []
    colors = {'summer':{'25':'#DE2C1B','24':'#D6806D','26':'#820707'}, 
              'winter' :{'25':'#2732CC','24':'#8295E0','26': "#200F80"},
              'spring' :{'25':"#00AE2E",'24':'#85F57F','26':'#046E35',},
              'autumn' :{'25':"#FFE202",'24':"#F5E6A1",'26':"#ADAA0C"}}
    for date in data.index:
        szn, yr = date.split('/')
        color = colors[szn][yr]
        clr_out.append(color)
        y.append(data[date].to_numpy())
        x.append(data.columns.to_numpy().astype(float))
        # coeff = np.polyfit(data.index.to_numpy().astype(float), data[date].to_numpy(), 5)
        # p = np.poly1d(coeff)
        # fit = (p(data.index.to_numpy().astype(float)))
        # fits.append(fit)
        leg.append(f'number distribution for {date}')# | {m:.2f}x + {b:.2f} | R2= {r:.4f} | corr = {cor:.2f}')
    scat_plot_Fact(x,y,leg,clr_out, append=append)


def scat_plot_dist(x,y,legs,colors,append = 0):
    plt.ion()
    fig, ax = plt.subplots()
    lines = []
    points = []
    for i in range(len(y)):
        L = ax.plot(x[i], y[i],label = legs[i], color = colors[i],)
        lines.append(L)
        pnt = [ps]
        pnt.extend(list(caps))
        pnt.append(list(bar)[0])
        points.append(pnt)
    leg = ax.legend()
    lined = dict()
    for legline, origline,point in zip(leg.get_lines(), lines, points):       
        legline.set_picker(10)  # 10 pts tolerance
        lined[legline] = origline[0], point
        # input(lined[legline])

    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline, point = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        for p in point:
            p.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_ylabel('Diameter [nm]')
    ax.set_xlabel('ss [%]')
    if append ==0:
        ax.set_title(f"Activation Diameter of CCN")
    else:
        ax.set_title(f"Activation Diameter of CCN for {append}")
    input('Press enter to exit plot...')
    plt.ioff()

if __name__ == '__main__':
    bad_dates = [[pd.to_datetime('10/01/2025 00:00:00'),pd.to_datetime('12/15/2025 00:00:00')],[pd.to_datetime('07/01/2025 00:00:00'),pd.to_datetime('07/24/2025 00:00:00')],[pd.to_datetime('08/06/2025 00:00:00'),pd.to_datetime('08/14/2025 00:00:00')]]
    cpc = r"C:\Users\bensy\Documents\Research\app_CPC.csv"
    smps =[r"C:\Users\bensy\Documents\Research\2024_SMPS_NumberSizeDist_1hr.csv",r"C:\Users\bensy\Documents\Research\SMPS_NumberSizeDist_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\2026_SMPS_NumberSizeDist_1hr.csv"]  #list(input('Provide paths to SMPS file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    ccn = [r"C:\Users\bensy\Documents\Research\CCN_Processed_2024_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2025_1hr.csv",r"C:\Users\bensy\Documents\Research\CCN_Processed_2026_1hr.csv"]     #list(input('Provide paths to CCN file(s). Seperate multiples with a comma: ').replace('"','').split(','))
    master =  r"C:\Users\bensy\Downloads\MasterDataFile_ChemAOPsCCNSMPSMET_June2024-Oct2025.csv"
    f = 'ME'
    dataout = r"C:\Users\bensy\Documents\Research\CCN_activation_diameter_test.csv"
    data,ss_cols,diam_cols = comb_files(smps,ccn,cpc, freq=f)
    smps_data = data[diam_cols]
    # Dist_curve(smps_data)
    data.to_csv(r"C:\Users\bensy\Documents\Research\check.csv")
    mask = pd.Series(False, index=data.index)
    for date in bad_dates:
        mask |= (data.index >= date[0]) & (data.index <= date[-1])
    data = data[~mask]
    codf = pd.DataFrame()
    actdf = pd.DataFrame()
    for row in range(len(data.index.to_numpy())):
        d = data.iloc[row]
        # cut_offs = find_cutoff(d,diam_cols,ss_cols)
        act_perc =find_activation(d,diam_cols,ss_cols)
        if row == 0:
            dt = pd.to_datetime(data.index.to_numpy()[row])
            year = str(dt.year).replace('20','')
            month = dt.month
            day = dt.day
            # codf = pd.DataFrame(cut_offs).T
            actdf = pd.DataFrame(act_perc)
            actdf = actdf.set_index('Diameter')
            if f == 'ME':
                # codf.columns = [f'D_cutoff_{month}/{year}',f'D_lower_{month}/{year}',f'D_upper_{month}/{year}']
                actdf = actdf.add_suffix(f'_{month}/{year}')
            else: 
                # codf.columns = [f'D_cutoff_{day}/{month}/{year}',f'D_lower_{day}/{month}/{year}',f'D_upper_{day}/{month}/{year}']
                actdf = actdf.add_suffix(f'_{day}/{month}/{year}')
        else: 
            dt = pd.to_datetime(data.index.to_numpy()[row])
            year = str(dt.year).replace('20','')
            month = dt.month
            day = dt.day
            # df = pd.DataFrame(cut_offs).T
            adf = pd.DataFrame(act_perc)
            adf = adf.set_index('Diameter')
            if f == 'ME':
                # df.columns = [f'D_cutoff_{month}/{year}',f'D_lower_{month}/{year}',f'D_upper_{month}/{year}']
                adf = adf.add_suffix(f'_{month}/{year}')
            else: 
                # df.columns = [f'D_cutoff_{day}/{month}/{year}',f'D_lower_{day}/{month}/{year}',f'D_upper_{day}/{month}/{year}']
                adf = adf.add_suffix(f'_{day}/{month}/{year}')
            # codf = pd.merge(codf,df, left_index=True, right_index=True)
            actdf = pd.merge(actdf, adf, left_index=True, right_index=True)
    # input(actdf)
    Fact_curve(actdf, freq='S')
    # cut_off_curve(codf, freq ='S')
    out = r"C:\Users\bensy\Documents\Research\CCN_activation_diameter.csv"#input("Enter filepath to export data as a csv, or press 'enter' to skip: ")
    if out != '':
        codf.to_csv(out)

