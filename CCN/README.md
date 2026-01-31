## Folder for CCN processing and generation of NASA AMES formated files along with general CSV files

### 1: CCN_main.py
  Main hub for CCN Processing. Run this script to process CCN data. Either explicity pass file paths prior to running or set "dev" variable to False to allow for inputs of file names to be passed during runtime.
### 2: CCN_process.py
  Functions called by CCN_main.py to process the data. Shouldn't need to be adjusted. All functions can be run independantly for smaller indepth data options
### 3: CCN_EBAS_convert.py
  Functions called to convert the CSV into a NASA AMES file in accordance with EBAS GWA formating
