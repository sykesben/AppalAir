"""
Date: 7/10/2026
Author: Ben Sykes
Purpose: Diurnal Trend Analyis
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from plotgen import monthly_box_call, hourly_box_call,line_plot, scat_call
from scipy.special import erf
pd.set_option('mode.chained_assignment', None)
# plt.rcParams['font.size'] = 18
# plt.rcParams['axes.titlesize'] =18

def smps_data_corr(files, freq='D'):
    """
    Clean and output the PSD produced by the SMPS

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    """
    smps_out = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        smps = pd.read_csv(f) #read in smps file\
        print(smps.columns)
        try:
            smps = smps.set_index("DateTime Sample Start") #Set index
        except:
            smps = smps.set_index("Datetime(UTC)") #Set index
        smps.index = pd.to_datetime(smps.index, format='mixed')

        numsmps = [s for s in smps.columns.to_numpy() if (('.' in s) and (s.split('.')[0].isdigit()))|(s.isdigit())]
        # numerical sort
        numsmps = sorted(numsmps, key=lambda x: float(x))
        smps = smps[numsmps]
        new_nums = {}
        for num in numsmps:
            try:
                # ensure all bin diameters are 2 decimal places to correctly merge
                float_num = float(num)
                new_nums[num] = f"{float_num:.2f}"
            except:
                new_nums[num] = num
        smps = smps.rename(columns=new_nums)
        #Remove the duplicate headers
        smps = smps.T.groupby(level=0).mean().T
        smps = smps.resample(freq).mean() 
        smps.index.names = ['Date']
        smps_out = smps if smps_out.empty else pd.concat([smps_out, smps])
    numsmps = np.asarray(list(new_nums.values()))
    smps_out = smps_out[numsmps]
    return smps_out, numsmps

def smps_means(files,freq='d'):
    '''
    Takes in a list of smps files and returns some values
    ----------

    Parameters
    ++++++++++
    files : [list of str] Paths to SMPS files
    freq : [str] Resample frequency for DataFrame

    Returns
    ++++++++++
    smps : [DataFrame] Combined SMPS data from all inputted files
    cols : [list of str] Names of used columns from SMPS output
    '''
    smps = pd.DataFrame()
    for i in range(len(files)): #read in smps files and combine
        f = files[i]
        file =pd.read_csv(f) #read in smps file
        try: 
            file = file.set_index('Datetime(UTC)')
        except:
            file=file.set_index("DateTime Sample Start") #Set index
        smps = file if smps.empty else pd.concat([smps, file])
    smps.index = pd.to_datetime(smps.index, format='mixed')
    smps = smps[smps['Total Concentration (#/cm³)'].notna()]
    cols = ['Median (nm)',"Mean (nm)",'Geo. Mean (nm)','Mode (nm)','Geo. Std. Dev','Total Concentration (#/cm³)']
    smps = smps[cols]
    return smps

def timeOfDay(hour):
    dscrpt= []
    if (hour<7)|(hour>19): # 7pm->7am
        dscrpt.append('night')
    else:
        dscrpt.append('day')
    if (hour<=8)&(hour>=6):
        dscrpt.append('dawn')
    if (hour<=19)&(hour>=17):
        dscrpt.append('dusk')
    if ((hour>=16)&(hour<=18))|(hour ==12):
        dscrpt.append('traffic')
    return dscrpt

def hourlyAvg(data):
    hours = data.index.hour.to_numpy(dtype = str)
    hours = np.asarray([('0'+m)[-2:] for m in hours])
    data['hour'] = hours
    hourly = data.groupby(['hour']).mean()
    return(hourly)

def timeSlct(data,times):
    mask = (data.index >= times[0]) & (data.index <= times[-1])
    out = data.copy().loc[mask]
    return out

def normData(data):
    norm = (data-data.min())/(data.max()-data.min())
    norm = norm.dropna(how='all')
    return norm

def standData(data):
    stand = (data-data.mean())/(data.std())
    stand = stand.dropna(how='all')
    return stand

SMPS_files = [r"C:\Users\bensy\Documents\Research\SMPS_2024\2024_SMPS_NumberSizeDist_1hr_clean.csv",
              r"C:\Users\bensy\Documents\Research\SMPS_2025\2025_SMPS_NumberSizeDist_1hr_clean.csv",
              r"C:\Users\bensy\Documents\Research\SMPS_2026\2026_SMPS_NumberSizeDist_1hr_clean.csv"]
file_out = r"C:\Users\bensy\Documents\Research\DayTime_SMPS.csv"
kappa_in = r"C:\Users\bensy\Documents\Research\CCN_SMPS_processed_parameters.csv"
march26 = [pd.to_datetime('03/01/2026'),pd.to_datetime('04/01/2026')]
june26 = [pd.to_datetime('06/01/2026'),pd.to_datetime('07/01/2026')]
spring26 = [pd.to_datetime('03/01/2026'),pd.to_datetime('06/01/2026')]
spring25 = [pd.to_datetime('03/01/2025'),pd.to_datetime('06/01/2025')]
summer25 = [pd.to_datetime('06/01/2025'),pd.to_datetime('09/01/2025')]
january25 = [pd.to_datetime('01/01/2025'),pd.to_datetime('02/01/2025')]
january26 = [pd.to_datetime('01/01/2026'),pd.to_datetime('02/01/2026')]
file = smps_means(SMPS_files)
kfile =pd.read_csv(kappa_in) #read in smps file
kfile = kfile.set_index('Datetime')
kfile.index = pd.to_datetime(kfile.index, format='mixed')
Dclow =kfile['Dcrit(ss=0.1)'].copy()
Dchigh =kfile['Dcrit(ss=0.7)'].copy()
geoMean = file['Geo. Mean (nm)'].copy()
file.index = pd.to_datetime(file.index) - pd.Timedelta('5h')
file.index.rename('Datetime(EST)',inplace=True)
file['hour'] = file.index.hour
freq_map = 'W'
PSD,bins = smps_data_corr(SMPS_files,freq=freq_map)
# input(PSD.T)
# largest = 
geoMean = geoMean.resample(freq_map).mean()
Dclow = Dclow.resample(freq_map).median()
Dchigh = Dchigh.resample(freq_map).median()
# input(f'{len(PSD.index.to_numpy())} vs {len(geoMean.index.to_numpy())}')
geoMean = geoMean.loc[np.isin(geoMean.index.date, PSD.index.date)]
Dclow = Dclow.loc[np.isin(Dclow.index.date, PSD.index.date)]
Dchigh = Dchigh.loc[np.isin(Dchigh.index.date, PSD.index.date)]
# input(geoMean)
bins= bins.astype(float)
PSD = PSD.clip(upper=5000)
PSD = PSD.fillna(0)

PSD.dropna(how='all',inplace=True)
dates = PSD.index.to_list()
PSD_data = PSD.T.to_numpy()
max_i = []
for col in PSD_data.T:
    if np.max(col)<15:
        max_i.append(max_i[-1])
    else:
        pull = np.where(col<100, bins, np.nan)
        pull = pull[~np.isnan(pull)]
        pull = np.where(pull>100,pull,np.nan)
        pull = pull[~np.isnan(pull)]
        try:
            first = pull[0]
        except:
            first = np.nan
        max_i.append(first)
fig, ax = plt.subplots()
# hm = ax.imshow(PSD_data, cmap='GnBu', interpolation='nearest')

# Add colorbar
mean = ax.plot(geoMean,color= '#310045',linestyle ='dashed',label =r'D$_{GeoM}$')
# ax.fill_between(Dclow.index.to_numpy(),Dclow.to_numpy(),Dchigh.to_numpy(),color= "#350A91", alpha =0.5)
dcl = ax.plot(Dclow,color= "#EFA105",linewidth=3,label =r'D$_{crit}$(ss=0.1)')
dch = ax.plot(Dchigh,color= "#C14A29",linewidth=3,label =r'D$_{crit}$(ss=0.7)')
# top = ax.plot(dates,max_i,color= '#E0634A',linewidth =5,alpha = 0.25)

hm = ax.pcolormesh(dates, bins, PSD_data,
                    cmap='YlGnBu',
                    shading='cool',
                    fc = 'red')
xt=dates[5]
yt=max_i[5]
# ax.annotate(r'N$_{CN}$<100',xy=(xt, yt), xycoords='data',xytext=(1.5, 30.5),color ='#E0634A',size=14,
#              textcoords='offset points',bbox=dict(boxstyle='round',facecolor="#E7CCC7",lw =0,alpha=0.25))
# xdl = Dclow.index.to_numpy()[2]
# ydl = Dclow.to_numpy()[2]
# ax.annotate(r'D$_{Crit}$(ss=0.15)',xy=(xdl, ydl), xycoords='data',xytext=(-75, 30.5),color ="#671985",size=14,
#              textcoords='offset points',bbox=dict(boxstyle='round',facecolor="#DDADF0",lw =0,alpha=0.25))
# xdh = Dchigh.index.to_numpy()[2]
# ydh = Dchigh.to_numpy()[2]
# ax.annotate(r'D$_{Crit}$(ss=0.7)',xy=(xdh, ydh), xycoords='data',xytext=(-75, 30.5),color ="#271257",size=14,
#              textcoords='offset points',bbox=dict(boxstyle='round',facecolor="#896DC6",lw =0,alpha=0.25))

# # xg = geoMean.index.to_numpy()[-25]
# yg = geoMean.to_numpy()[-25]
# ax.annotate(r'D$_{GeoM}$',xy=(xg, yg), xycoords='data',xytext=(1.5, -40.5),color = '#310045',size=15,
#              textcoords='offset points',bbox=dict(boxstyle='round',facecolor="#DCCBE3",lw =0,alpha=0.25))
ax.set_yscale('log')
cbar = ax.figure.colorbar(hm, ax=ax)
cbar.ax.set_ylabel('Particle Concentration', rotation=-90, va="bottom")
# Major tick every 3 months
ax.yaxis.set_label_text('Diameter [log(nm)]')
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax.legend()

# Format as YYYY-MM
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))

plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
# # Show ~10 x-axis labels
# xticks = np.linspace(0, len(dates)-1, 10, dtype=int)
# ax.set_xticks(xticks)
# ax.set_xticklabels(np.array(dates)[xticks], rotation=45, ha='right')

# # Show ~10 y-axis labels
# yticks = np.linspace(0, len(bins)-1, 10, dtype=int)
# ax.set_yticks(yticks)
# ax.set_yticklabels(np.array(bins)[yticks])
# dates = PSD.index.to_list()

plt.title("Change in Size Distribution")
plt.show()
input()


mar26 = timeSlct(file,march26)
jun26 = timeSlct(file,june26)
spr26 = timeSlct(file,spring26)
spr25 = timeSlct(file,spring25)
sum25 = timeSlct(file,summer25)
jan26 = timeSlct(file,january26)
jan25 = timeSlct(file,january25)

ToD = file.hour.apply(timeOfDay).to_numpy()
diurnal = []
traffic = []
DorD =[]
for i in ToD:
    diurnal.append(i[0])
    if 'traffic' in i:
        traffic.append('traffic')
    else:
        traffic.append('False')
    if 'dusk' in i:
        DorD.append('dusk')
    elif 'dawn' in i:
        DorD.append('dawn')
    else:
        DorD.append("False")
file['dayNight'] = diurnal
file['dawnDusk'] = DorD
file['traffic'] = traffic

dayfile = file.copy()[file['dayNight'] == 'day']
nightfile = file.copy()[file['dayNight'] == 'night']
trafficFile = file.copy()[file['traffic'] == 'traffic']
print(len(file))
print(len(dayfile))
print(len(nightfile))

hmar26 = hourlyAvg(mar26)
hjun26 = hourlyAvg(jun26)
hspr26 = hourlyAvg(spr26)
hspr25 = hourlyAvg(spr25)
hsum25 = hourlyAvg(sum25)
hjan25 = hourlyAvg(jan25)
hjan26 = hourlyAvg(jan26)
hspr26.to_csv(r"C:\Users\bensy\Documents\Research\Spring_diurnal_trends_SMPS.csv")
nmar26 = standData(hmar26)
njun26 = standData(hjun26)
nspr26 = standData(hspr26)
nspr25 = standData(hspr25)
nsum25 = standData(hsum25)
njan26 = standData(hjan26)
njan25= standData(hjan25)

dayfile.to_csv(file_out)
# input(dayfile.columns.to_numpy())
monthly_box_call(file,'Total Concentration (#/cm³)',y_label=r'N$_{CN}$[#/cm3]',title = 'Monthly SMPS Concentration')
monthly_box_call(file,'Geo. Mean (nm)',y_label=r'Geo.Mean D(nm)',title = 'Monthly Geo. Mean Diameter')
# monthly_box_call([file,trafficFile], 'Kappa(ss=0.25)',y_label='Κ(ss=0.25)',title = 'Traffic Monthly Hygroscopicity Data',keys=['Total','Traffic'])
# hourly_box_call(spr25, 'Geo. Mean (nm)',y_label='Geo. Mean D(nm)',title = 'Spring 26 Hourly Hygroscopicity Data')
y = [nspr26['Kappa(ss=0.25)'].to_numpy(),nspr25['Kappa(ss=0.25)'].to_numpy(),nsum25['Kappa(ss=0.25)'].to_numpy()]#nspr26["org/total"].to_numpy(),nspr26['SO4/total'].to_numpy(),nspr26['NO3/total'].to_numpy(),]
org = [nspr26['org/total'].to_numpy(),nspr25['org/total'].to_numpy(),nsum25['org/total'].to_numpy(),
     njan25['org/total'].to_numpy(),njan26['org/total'].to_numpy()]
OOA = [nspr26['NO3/total'].to_numpy(),nspr25['NO3/total'].to_numpy(),nsum25['NO3/total'].to_numpy()]
tot = [nspr26['Total Concentration (#/cm³)'].to_numpy(),nspr25['Total Concentration (#/cm³)'].to_numpy(),nsum25['Total Concentration (#/cm³)'].to_numpy()]
martwnsx= [nmar26['Kappa(ss=0.25)'].to_numpy(),nmar26['Total mass'].to_numpy(),nmar26['SO4[ug/m3]'].to_numpy(),nmar26['Org[ug/m3]'].to_numpy()]
aprtwnsx= [njun26['Kappa(ss=0.25)'].to_numpy(),njun26['Total Concentration (#/cm³)'].to_numpy(),njun26['SO4[ug/m3]'].to_numpy(),njun26['Org[ug/m3]'].to_numpy()]
twntysx= [nspr26['Kappa(ss=0.25)'].to_numpy(),nspr26['Total mass'].to_numpy(),nspr26['SO4[ug/m3]'].to_numpy(),nspr26['Org[ug/m3]'].to_numpy()]
twntyfv= [nspr25['Kappa(ss=0.25)'].to_numpy(),nspr25['Total Concentration (#/cm³)'].to_numpy(),nspr25['SO4[ug/m3]'].to_numpy(),nspr25['Org[ug/m3]'].to_numpy()]
twntysum= [nsum25['Kappa(ss=0.25)'].to_numpy(),nsum25['Total Concentration (#/cm³)'].to_numpy(),nsum25['SO4[ug/m3]'].to_numpy(),nsum25['Org[ug/m3]'].to_numpy()]
jantwsx= [njan26['Kappa(ss=0.25)'].to_numpy(),njan26['Total Concentration (#/cm³)'].to_numpy(),njan26['SO4[ug/m3]'].to_numpy(),njan26['Org[ug/m3]'].to_numpy()]
legTot = ['Κ(ss=0.25)',r'Mass',r'${SO4}$',r'${org}$']
x = hspr25.index.to_numpy()
legs = ['Spring 26','Spring 25','Summer 25']#['K(ss=0.25)',r'F$_{OA}$',r'F$_{SO4}$',r'F$_{NO3}$']
# slct= {'Kappa(ss=0.25)':'SO4/total'}
# scat_call(hspr26,slct,x_label=r'F$_{SO4}$', y_label=r'Κ(ss=0.25)',title='Hygroscopicity vs Sulfate Fraction Spring 26 Hourly')
        # 's26 Geo.Mean D(nm)','s25 Geo.Mean D(nm)']
# line_plot(x,y,legs,x_label='hour of day',y_label='Κ(ss=0.25)[1]', title='Daily Trend of Seasons')
# line_plot(x,OOA,legs,x_label='hour of day',y_label=r'F$_{NO3}$[1]', title='Daily Trend of Seasons')
# line_plot(x,tot,legs,x_label='hour of day',y_label=r'N$_{CN}$[1]', title='Daily Trend of Seasons')
line_plot(x,martwnsx,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='March 26 Time of Day')
# line_plot(x,aprtwnsx,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='June 26 Time of Day')
line_plot(x,twntysx,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='Spring 26 Time of Day')
line_plot(x,twntyfv,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='Spring 25 Time of Day')
line_plot(x,twntysum,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='Summer 25 Time of Day')
line_plot(x,jantwsx,legTot,x_label='hour of day',y_label=r'normalized trends[1]', title='Jan 26 Time of Day')
# monthly_box_call(nightfile, 'Kappa(ss=0.25)',y_label='Κ(ss=0.25)',title = 'Nocturnal Monthly Hygroscopicity data')

    
input(dayfile)