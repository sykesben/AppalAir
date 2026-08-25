#Aydan Gibbs
#1/22/25
#box and wisker plots for smps


#Next Steps
#sorting files by thier names to put them in order on the graph

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import smps
import os
import requests
from pathlib import Path


def main(): #skipping i rows for correct data range (3948 for June 2025 onward)
    import pandas as pd

    # Load the CSV file into a DataFrame
    file_path = 'C:\\Users\\ze_ba\\OneDrive\\Desktop\\AppalAIR\SMPS-Data\\testing\whisker\\SMPS_NumberSizeDist_2025_1hr.csv'  # Replace with your CSV file path
    df = pd.read_csv(file_path, skiprows = 3948)

    # Define the range of columns to sum
    # Assuming your columns are named 'col_1', 'col_2', ..., 'col_500'
    start_col_index = 160  # Index of the first column (9 for 13.47) (121 for 100nm) (160 for 200nm)
    end_col_index = 236  # Index of the last column (e.g., 'col_500')

    # Select the columns within the specified range
    # Using iloc for integer-location based indexing
    columns_to_sum = df.iloc[:, start_col_index : end_col_index + 1]

    # Sum the values across the selected columns for each row
    # axis=1 specifies summing across columns (row-wise)
    df['summed_columns'] = columns_to_sum.sum(axis=1)
    df['DateTime Sample Start'] = pd.to_datetime(df.iloc[:, 0])
    # Display the DataFrame with the new summed column
    print(df['summed_columns'])
    df.plot(x = 'DateTime Sample Start', y = 'summed_columns', title='Number of particles from 100nm to 800nm from June 2025 through October 2025',
            xlabel='Time', ylabel='#',
            color='skyblue')
    plt.show()

if __name__:
    main()