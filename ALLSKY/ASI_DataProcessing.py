import os
import pandas as pd
import numpy as np
import sys
from pandas.tseries.offsets import DateOffset
import datetime

pd.options.mode.chained_assignment = None  # removes a warning that this code generates that doesn't matter

np.set_printoptions(threshold=sys.maxsize)

os.chdir('C:/Users/lidar/Desktop/FCT22 Data/DATA/FCT22 Data/RESULTS/FILES/asi_16585')

try:
	year = int(input("Enter the current year as an integer: ")) #The year of the month you are averaging
	endMonth = int(input("Enter the month of the latest csv as an integer: ")) #The month you are averaging/the month of the latest asi csv
	endDay = int(input("Enter the day of the most recent csv as an integer: ")) #The day of the latest asi csv
except Exception as e:
        print(e)


td = datetime.date(year,1,1) #start timeframe at beginning of year
tf = datetime.date(year,endMonth,endDay) #end timeframe on most recent csv
delta = (tf-td).days + 1 #+1 to include the last day
print("Numbr of days: " + str(delta))
df_ASI = pd.DataFrame()

#Write in every csv to a df
#Remove extra headers
#separate the first column by semicolon if need be
#remove columns besides date, time and cloudiness
#combine date and time into one pandas datetime column
#combine all data from each year into 1 csv
#Save Unaveraged and Hourly Averaged data to a csv per year

for day in range(delta): #loop through every day in the timeframe
    #turn current day into the proper filename format for reading in the csv
    date = (td + datetime.timedelta(days = day)).isoformat()
    date = date.replace("-", "")
    file = "asi_16585_" + date + ".csv"
    
    #try to read in the csv for the day
    try:
        df = pd.read_csv(file)
    except Exception as e:
        print(e)
    
    #Cleans up the csv if it came out of the program with all the data in one column
    if (df.columns[0].find(";") > 0): #If the first column header contains a semicolon split the csv by semicolon
        columnnames = list(df.columns)[0].split(';')
        dfsplit = df.iloc[:,0].str.split(";", expand=True)
        dfsplit.columns = columnnames
        dfnext = dfsplit
    else:
        dfnext = df
      
    #removes any extra headers in the csvs (if there is a "time" value in the time column)
    rows = dfnext[dfnext["time"] == "time"].index 
    dfnext.drop(rows, inplace=True)
    
    dfclean = dfnext[['%date', 'time', 'Cloudiness']]

    try:
        dfclean = dfclean.astype({'%date' : str, 'time' : str, 'Cloudiness' : float}) #convert date and time column to string for next line of code
        dfclean['Datetime'] = pd.to_datetime(dfclean['%date'] + ' ' + dfclean['time'], errors='coerce') #create the DateTime column by combining the date and time columns
        dfout = dfclean.drop(['%date', 'time'], axis=1) #drop the date and time columns
        dfout = dfout[['Datetime', 'Cloudiness']] #reorder the dataframe so the DateTime column goes first
        dfout = dfout.set_index('Datetime') #set the datetime column to the index, gets rid of the numbered indeces
    except Exception as e:
        print(f"Error combining date and time: {e}")
        
    df_ASI = pd.concat((df_ASI, dfout))

#Write out unaveraged csv
df_ASI.to_csv("UNAVG_ASI_"+ str(year) +"_CloudFraction.csv")
print("UNAVG_ASI_"+ str(year) +"_CloudFraction.csv outputted.")
#Write out averaged csv
df_ASI.resample('h').mean().to_csv("HRAVG_ASI_"+ str(year) +"_CloudFraction.csv")
print("HRAVG_ASI_"+ str(year) +"_CloudFraction.csv outputted.")

input("Press Enter to Exit.")