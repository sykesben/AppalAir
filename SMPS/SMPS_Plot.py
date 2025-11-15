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

def main():

    dataTotal = pd.DataFrame()
    dataType = "Number"
    fileNames = []
    csvpath = Path(input("\nInput the full path of the csv youd like to access:\n"))
    print(csvpath)
    df = pd.read_csv(csvpath)


    columns_to_sum = df.iloc[:, 11 : 263 + 1]

    df['sum_of_100_nm_onward'] = df[columns_to_sum].sum(axis=1)
    print(df.head())


def fuckingAI():
    import pandas as pd

    # Load the CSV file into a DataFrame
    file_path = 'C:\\Users\\ze_ba\\OneDrive\\Desktop\\AppalAIR\SMPS-Data\\testing\whisker\\SMPS_NumberSizeDist_2025_1hr.csv'  # Replace with your CSV file path
    df = pd.read_csv(file_path)

    # Define the range of columns to sum
    # Assuming your columns are named 'col_1', 'col_2', ..., 'col_500'
    start_col_index = 11  # Index of the first column (e.g., 'col_1')
    end_col_index = 216  # Index of the last column (e.g., 'col_500')

    # Select the columns within the specified range
    # Using iloc for integer-location based indexing
    columns_to_sum = df.iloc[:, start_col_index : end_col_index + 1]

    # Sum the values across the selected columns for each row
    # axis=1 specifies summing across columns (row-wise)
    df['summed_columns'] = columns_to_sum.sum(axis=1)

    # Display the DataFrame with the new summed column
    print(df['summed_columns'])

    # You can also save the result to a new CSV file
    # df.to_csv('output_with_summed_columns.csv', index=False)

if __name__:
    fuckingAI()