This Code requires a folder that you will have to put the csv’s from the FCT software into, copy this folder path into the os.chdir(''); line of code

Takes 3 inputs: The current year (int), the month of the most recent csv (int), the day of the most recent csv(integer)

It will output an hourly averaged csv and an unaveraged csv for the year. The ASI does not take data overnight so there will be gaps in the hourly averaged data during these times. The unaveraged data should have no gaps but will not have consistent time deltas between each row.

