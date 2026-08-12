"""
Date: 7/10/2026
Author: Ben Sykes
Purpose: Diurnal Trend Analyis
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from CombineData import comb_files,smps_means
from plotgen import monthly_box_call, hourly_box_call,line_plot, scat_call,ss_box_call
from scipy.special import erf
pd.set_option('mode.chained_assignment', None)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.titlesize'] =18

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

def ssMelt(data, kept):
    cols = [col for col in data.columns.to_numpy() if kept in col]
    data = data[cols]
    df_long = (data.reset_index(drop=True).melt( value_name=kept,var_name='ss'))
    # Extract the supersaturation value
    df_long['ss'] = (df_long['ss'].str.extract(r'Kappa\(ss=([\d.]+)\)').astype(float))
    df_long = df_long.set_index('ss')
    return df_long

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

file_in = r"C:\Users\bensy\Documents\Research\CCN_SMPS_chem.csv"
file_out = r"C:\Users\bensy\Documents\Research\DayTime_CCN_SMPS_processed_parameters.csv"
march26 = [pd.to_datetime('03/01/2026'),pd.to_datetime('04/01/2026')]
june26 = [pd.to_datetime('06/01/2026'),pd.to_datetime('07/01/2026')]
spring26 = [pd.to_datetime('03/01/2026'),pd.to_datetime('06/01/2026')]
spring25 = [pd.to_datetime('03/01/2025'),pd.to_datetime('06/01/2025')]
summer25 = [pd.to_datetime('06/01/2025'),pd.to_datetime('09/01/2025')]
january25 = [pd.to_datetime('01/01/2025'),pd.to_datetime('02/01/2025')]
january26 = [pd.to_datetime('01/01/2026'),pd.to_datetime('02/01/2026')]
file =pd.read_csv(file_in)
file = file.set_index('Datetime')
file.index = pd.to_datetime(file.index) - pd.Timedelta('5h')
file.index.rename('Datetime(EST)',inplace=True)
file['hour'] = file.index.hour

mar26 = timeSlct(file,march26)
jun26 = timeSlct(file,june26)
spr26 = timeSlct(file,spring26)
spr25 = timeSlct(file,spring25)
sum25 = timeSlct(file,summer25)
jan26 = timeSlct(file,january26)
jan25 = timeSlct(file,january25)

kspr25 = ssMelt(sum25,'Kappa')
ss_box_call(kspr25,'Kappa',y_label=r'k$_{CCN}$',title = 'Spring kappa curve')

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
hspr26.to_csv(r"C:\Users\bensy\Documents\Research\Spring_diurnal_trends_2.csv")
nmar26 = standData(hmar26)
njun26 = standData(hjun26)
nspr26 = standData(hspr26)
nspr25 = standData(hspr25)
nsum25 = standData(hsum25)
njan26 = standData(hjan26)
njan25= standData(hjan25)

dayfile.to_csv(file_out)
# input(dayfile.columns.to_numpy())

monthly_box_call(file,'Kappa(ss=0.1)',y_label=r'k$_{CCN}$(ss=0.1)',title = r'Monthly k$_{CCN}$')
monthly_box_call(file,'Kappa(ss=0.25)',y_label=r'k$_{CCN}$(ss=0.25)',title = r'Monthly k$_{CCN}$')
monthly_box_call(file,'Kappa(ss=0.7)',y_label=r'k$_{CCN}$(ss=0.7)',title = r'Monthly k$_{CCN}$')
monthly_box_call(file,'Geo. Mean (nm)',y_label=r'Geo.Mean D(nm)',title = 'Monthly Geo. Mean Diameter')
# monthly_box_call(file,'N(cm-3)_cor_setpt0.25',y_label=r'N$_{CCN}$[#/cm3](ss=0.25)',title = 'Monthly CCN Concentration Data')
# monthly_box_call([dayfile,nightfile],'N(cm-3)_cor_setpt0.25',y_label=r'N$_{CCN}$[#/cm3](ss=0.25)',title = 'Monthly CCN Concentration Data',keys=['day','night'])
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