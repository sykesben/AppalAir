This Code requires a folder that you will have to put the csv’s from the FCT software into, copy this folder path into the os.chdir(''); line of code

You will have to change the end date in the code to line up with the date of the most recent csv from the FCT software. The start date should always be the first day of the year that you are in

It will output an hourly averaged csv and an unaveraged csv for the year. The ASI does not take data overnight so there will be gaps in the hourly averaged data during these times. The unaveraged data should have no gaps but will not have consistent time deltas between each row.

