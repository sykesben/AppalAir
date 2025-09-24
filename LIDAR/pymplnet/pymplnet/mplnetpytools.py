# mplnetpytools.py
# Collection of tools for working with MPLNET data
# Author: Jordan Greene

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import requests
import bs4
import re
import os
import netCDF4 as nc
from mplnetaccess import MPLNetAccess
import warnings

# TODO Build into class to hold user/pass for needed access
def get_mplnet_html(url):
    """
    Returns the html for the mplnet address given
    Used to scrape for file variables to assist in downloading
    https://mplnet.gsfc.nasa.gov/out/data/V3_partners/Appalachian_State/
    
    Args:
        url (str): The url for the mplnet html
    
    Returns:
        requests.Response object: The html for the mplnet address given
    """

    # Get username and password from encrypted file
    key = MPLNetAccess()
    key.getAccess()

    # creates a requests session
    s = requests.Session()

    # authenticates the session
    s.auth = (key.getUser(), key.getPass())

    # Consider just returning the text
    # return s.get(url).text

    # returns the html object for the data
    return s.get(url)

def get_mplnet_data(url, file_path_name):
    """
    Writes the mplnet data from file at url to file_path_name

    Args:
        url (str): The url for the mplnet data file
        file_path_name (str): The full path and name of the file to be saved

    Returns:
        Boolean value indicating if the file was downloaded and saved successfully
    """

    # Download the file with request

    # Get username and password from encrypted file
    key = MPLNetAccess()
    key.getAccess()

    # creates a requests session
    s = requests.Session()

    # authenticates the session
    s.auth = (key.getUser(), key.getPass())
    response = s.get(url, stream=True)

    if(response.status_code != 200):
        return False

    # If file already exists, return True, don't download again
    if (os.path.exists(file_path_name)):
        return True
    else: 
        with open(file_path_name, 'wb') as f:
            for chunk in response.iter_content(chunk_size = 10 * 1024):
                if chunk:
                    f.write(chunk)

    # Verify file was downloaded
    return os.path.exists(file_path_name)

def parse_year(html):
    """
    Returns a list of the years available for a given site

    Args:
        html (str): The html for the mplnet address given

    Returns:
        list: The years available for a given site
    """

    # parses the html
    soup = bs4.BeautifulSoup(html, 'html.parser')

    # finds all the years available
    years = soup.find_all('a', {'href': re.compile(r'Y(\d){4}')})

    # returns the years from 2021 onwards as a list without the 'Y' or '/' from directory structure 
    years_strip = [year.text.strip('Y').strip('/') for year in years]
    return years_strip[5:]

def parse_month(html):
    """
    Returns a list of the months available for a given site and year

    Args:
        html (str): The html for the mplnet address given
    
    Returns:
        list: The months available for a given site and year
    """

    # parses the html
    soup = bs4.BeautifulSoup(html, 'html.parser')

    # finds all the months available
    months = soup.find_all('a', {'href': re.compile(r'M(\d){2}')})

    # returns the months as a list without the 'M' or '/' from directory structure
    return [month.text.strip('M').strip('/') for month in months]

def parse_day(html):
    """
    Returns a list of the days available for a given site, year, and month

    Args:
        html (str): The html for the mplnet address given
    
    Returns:
        list: The days available for a given site, year, and month
    """

    # parses the html
    soup = bs4.BeautifulSoup(html, 'html.parser')

    # finds all the days available
    days = soup.find_all('a', {'href': re.compile(r'D(\d){2}')})

    # returns the days as a list without the 'D' or '/' from directory structure
    return [day.text.strip('D').strip('/') for day in days]

def parse_file(html):
    """
    Returns a list of the files available for a given site, year, month, and day

    Args:
        html (str): The html for the mplnet address given
    
    Returns:
        list: The files available for a given site, year, month, and day
    """

    # parses the html
    soup = bs4.BeautifulSoup(html, 'html.parser')

    # finds all the files available
    files = soup.find_all('a', {'href': re.compile(r'MPL.+Appalachian_State\.nc4')})

    # returns the files as a list
    return [file.text for file in files]

def directory_select(name, options):
    """
    Returns user selected directory from a list

    Depreciated as GUI will guide user selection

    Args:
        name (str): The name of the directory
        options (list): The options for the directory

    Returns:
        str: The user selected directory
    """

    user_input = ''
    input_message = '\nChoose ' + name + ': '
    for option in options:
        input_message += option + ' '
    print()
    input_message += '\nYour choice: '
    while user_input.lower() not in options:
        user_input = input(input_message)

    print('You picked: ' + user_input)
    return user_input

def file_select(file_list):
    """
    Returns user selected file from a list
    
    Depreciated as GUI will guide user selection

    Args:
        file_list (list): The list of files
        
    Returns:
        str: The user selected file
    """

    user_input = ''
    input_message = '\nChoose file: \n'

    for index, file in enumerate(file_list):
        input_message += f'{index+1}) {file}\n'

    input_message += 'Your choice: '

    while user_input not in range(1, len(file_list)+1):
        user_input = int(input(input_message))

    print('You picked: ' + file_list[user_input-1])

    return file_list[user_input-1]

def create_directory(path):
    """
    Returns the url for the mplnet file path

    Args:
        path (str): The path to be folder to be created

    Returns:
        None
    """

    # If directory does not exist, create it
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        print('Creation of the directory %s failed' % path)

def get_user_path():
    """
    Returns the user specified path for the mplnet data

    Args:
        None

    Returns:
        str: The user specified local path for the mplnet data
    """

    file_path = input('Enter the path to the mplnet data: ')

    while (not os.path.exists(file_path)):
        print('No such directory')
        file_path = input('Enter the path to the mplnet data: ')

    return file_path

def leap_year(year):
    """
    Returns true if year is a leap year

    Args:
        year (int): The year to be checked

    Returns:
        bool: True if year is a leap year
    """

    return ((year % 4) or (year % 100 < 1) and (year % 400)) < 1

def buildSelectionList(years):
    """
    Returns a list of months, days, files, levels for selection

    Args:
        years (list): The years available for a given site

    Returns:
        months (list): The months available for a given site and year
        days (list): The days available for a given site, year, and month
        fileTypes (list): The files available for a given site, year, month, and day
    """

    # make sure years is not empty
    if not years:
        return []

    # Create list of months
    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    # Create list of days
    # Check if leap year and add extra day to February if true
    leapYear = int(max([leap_year(x) for x in list(map(int, years))]))
    # daysInMonths = [31, 28 + leapYear, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days = [str(x).zfill(2) for x in range(1, 32)]

    # Create list of file types
    fileType = ['NRB', 'CLD', 'AER', 'PBL']

    #create list of levels
    levels = ['L1', 'L15']

    return months, days, fileType, levels

# Class to build and store file variables
class FileVariables():
    """ 
    Class to build and store file variables
    
    Attributes:
        filetype (str): The file type
        filevars (list): The file variables
        nextIter (int): The next file variable to be filled
        selectedFileType (str): The selected file type
        selectedFileVars (list): The selected file variables
        
    Methods:
        setFileTypes(ft): Store the file types from the downloaded files
            into the filetype variable to allow user selection
        next(): Returns the next file variable to be filled
            based upon the user input
        peakNext(): Returns the next file variable to be filled
            based upon the user input
        storeCurrent(value): Returns the current file variable to be filled
            based upon the user input
        setFileVars(file): Returns the file variables for the selected file type
        printSelected(var): Prints the selected file variable
    """

    def __init__(self):
        # Initialize variables that hold file variables
        # Select file type from downloaded file types
        self.filetype = []
        self.filevars = {}
        self.nextIter = 0

        # Initialize variables that hold user selected file variables
        self.selectedFileType = ''
        # This will be netcdf variable
        self.selectedFileVars = ''

    def setFileTypes(self, ft):
        """
        Store the file types from the downloaded files
        into the filetype variable to allow user selection
        
        Args:
            ft (list): The list of file types
        
        Returns:
            None
        """
            
        self.filetype = ft

    def next(self):
        """
        Returns the next file variable to be filled 
        based upon the user input
        
        Args:
            None
        
        Returns:
            The next file variable to be filled
        """

        if self.nextIter == 0:
            self.nextIter += 1
            return self.filetype
        elif self.nextIter == 1:
            self.nextIter += 1
            return self.filevars
        else:
            return None

    def peakNext(self):
        """
        Returns the next file variable to be filled
        based upon the user input
        
        Args:
            None
        
        Returns:
            Boolean value indicating if there is a next file variable
        """

        if self.nextIter == 0:
            return True
        elif self.nextIter == 1:
            return True
        else:
            return False

    def storeCurrent(self, value):
        """
        Stores the current file variable to be filled 
        based upon the user input and the nextIter variable
        
        Args:
            None
        
        Returns:
            None
        """

        if self.nextIter == 1:
            self.selectedFileType = value
        elif self.nextIter == 2:
            self.selectedFileVars = value
        else:
            return None

    def setFileVars(self, file):
        """
        Opens the netcdf file and extracts the variables and populates
        the filevars dictionary
        
        Args:
            files (str): A full path to the netcdf file to be opened
        
        Returns:
            None
        """
        
        # Open the netcdf file and extract the variables into a dictionary for use
        with nc.Dataset(file, 'r') as ncFile:
            for name, variable in ncFile.variables.items():
                self.filevars[name] = {
                    'attributes': {attr_name: getattr(variable, attr_name) for attr_name in variable.ncattrs()},
                    'dimensions': variable.dimensions,
                    'shape': variable.shape,
                    'data': variable[:].copy()
                }


    def printSelected(self, var):
        """
        Prints the selected file variable

        Args:
            var (str): The selected file variable
            
        Returns:
            String with the selected file variable
        """

        return f'\nVariable: {var}\
                \n  Dimensions: {self.filevars[var]["dimensions"]}\
                \n  Shape: {self.filevars[var]["shape"]}'

# Class to build and store selection variables
class SelectionVariables():
    """
    Class to build and store selection variables

    Attributes:
        years (list): The list of years
        months (list): The list of months
        days (list): The list of days
        fileTypes (list): The list of file types
        level (list): The list of levels
        selectedYears (list): The list of selected years
        selectedMonths (list): The list of selected months
        selectedDays (list): The list of selected days
        selectedFileTypes (list): The list of selected file types
        selectedLevel (list): The list of selected level
        nextIter (int): The next selection variable to be filled

    Methods:
        getVars(html): Returns the years, months, days, and file types from the html
        next(): Returns the next selection variable to be filled
            based upon the user input
        peakNext(): Returns the next selection variable to be filled
            based upon the user input
        storeCurrent(value): Returns the current selection variable to be filled
            based upon the user input
        prepDownload(): Returns the urls, directories, and files to be downloaded
        download(url, dir, file): Downloads the selected files
        checkSelection(): Checks to make sure a selection is made in each category
        reset(): Resets the nextIter variable to 0
        printSelected(): Prints the selected variables

    TODO: Add in Level 1.0 and 1.5 options
    """

    def __init__(self):
        # Initialize variables that hold selection options
        self.years = [] 
        self.months =  []
        self.days = []
        self.fileTypes = []
        self.level = []
        self.nextIter = 0

        # Initialize variables that hold user selected options
        self.selectedYears = [] 
        self.selectedMonths =  []
        self.selectedDays = []
        self.selectedFileTypes = []
        self.selectedLevel = []

    def getVars(self, html):
        """
        Builds the years, months, days, and file types from the html
        to be used in the selection of files to be downloaded
        
        Args:
            html (str): The html from the mplnet website
            
        Returns:
            None    
        """
        
        self.years = parse_year(html)
        self.months, self.days, self.fileTypes, self.level = buildSelectionList(self.years)

    def next(self):
        """
        Returns the next selection variable to be filled 
        based upon the user input
        
        Args:
            None
        
        Returns:
            The next selection variable to be filled or None
        """

        if self.nextIter == 0:
            self.nextIter += 1
            return self.years
        elif self.nextIter == 1:
            self.nextIter += 1
            return self.months
        elif self.nextIter == 2:
            self.nextIter += 1
            return self.days
        elif self.nextIter == 3:
            self.nextIter += 1
            return self.fileTypes
        elif self.nextIter == 4:
            self.nextIter += 1
            return self.level
        else:
            return None

    def storeCurrent(self, value):
        """
        Returns the current selection variable to be filled 
        based upon the user input
        
        Args:
            value (list): The current selection values to be stored
        
        Returns:
            None
        """

        if self.nextIter == 1:
            self.selectedYears = value
        elif self.nextIter == 2:
            self.selectedMonths = value
        elif self.nextIter == 3:
            self.selectedDays = value
        elif self.nextIter == 4:
            self.selectedFileTypes = value
        elif self.nextIter == 5:
            self.selectedLevel = value
        else:
            return None

    def peakNext(self):
        """
        Returns the next selection variable to be filled
        based upon the user input
        
        Args:
            None
        
        Returns:
            Boolean value indicating if there is a next selection variable
        """

        if self.nextIter == 0:
            return True
        elif self.nextIter == 1:
            return True
        elif self.nextIter == 2:
            return True
        elif self.nextIter == 3:
            return True
        elif self.nextIter == 4:
            return True
        else:
            return False

    def reset(self):
        """
        Resets the nextIter variable to 0 and clears the selected values
        
        Args:
            None
        
        Returns:
            None
        """

        # reset selected values
        self.nextIter = 0
        self.selectedYears = []
        self.selectedMonths = []
        self.selectedDays = []
        self.selectedFileTypes = []
        self.selectedLevel = []

    def printSelected(self):
        """
        Prints the selected variables
        
        Args:
            None
        
        Returns:
            String with the selected variables
        """

        return 'Years selected: \n' + ' '.join(self.selectedYears)\
            + '\n\nMonths selected: \n' + ' '.join(self.selectedMonths)\
            + '\n\nDays selected: \n' + ' '.join(self.selectedDays)\
            + '\n\nFileTypes selected: \n' + ' '.join(self.selectedFileTypes)\
            + '\n\nLevel selected: \n' + ' '.join(self.selectedLevel) + '\n'

    # TODO: Remove prepDownload and download methods from this class
    def prepDownload(self):
        """
        Prepares a tuple the selected files to be downloaded including
        the urls, directories, and files
        
        Args:
            None
        
        Returns:
            Tuple of lists containing the urls, dirs, and, files to be downloaded
        """

        # Build list of files to be downloaded
        files = buildFileList(self.selectedYears, self.selectedMonths, self.selectedDays, self.selectedFileTypes, self.selectedLevel)

        # Build list of directories to be downloaded
        dirs = buildDirList(self.selectedYears, self.selectedMonths, self.selectedDays, self.selectedFileTypes, self.selectedLevel)

        # Build list of urls to be downloaded
        urls = buildURLList(self.selectedYears, self.selectedMonths, self.selectedDays, self.selectedFileTypes, self.selectedLevel)

        # Reset the selection variables
        return urls, dirs, files

    def download(self, url, dir, file):
        """
        Downloads the selected files
        
        Args:
            url (str): The url for the mplnet data file
            dir (str): The full path of the file to be saved, not including the file name
            file (str): The name of the file to be saved
        
        Returns:
            String with the file name and status of the download
        """
        
        # Download the files
        create_directory(dir)
        if os.path.exists(dir + file):
            return file + ' already exists'
        if get_mplnet_data(url, dir + file):
            return file + ' download successful'
        else: 
            return file + ' download failed'

    def checkSelection(self):
        """
        Checks to make sure a selection is made in each category
        
        Args:
            None
        
        Returns:
            Boolean value indicating if a selection is made in each category
        """

        return (self.selectedYears and self.selectedMonths and self.selectedDays and self.selectedFileTypes and self.selectedLevel)
            
def buildFileList(years, months, days, fileTypes, levels):
    """
    Returns a list of files to be downloaded based upon the user input
    
    Args:
        years (list): The list of years
        months (list): The list of months
        days (list): The list of days
        fileTypes (list): The list of file types
    
    Returns:
        files (list): The list of files to be downloaded
    """

    # Create list of files
    files = []
    for year in years:
        for month in months:
            for day in days:
                for fileType in fileTypes:
                    for level in levels:
                        files.append('MPLNET_V3_' + level + '_' + fileType + '_' + year + month + day + '_MPL44201' + '_Appalachian_State.nc4')

    return files

def buildDirList(years, months, days, fileType, levels, data_path = '../data/'):
    """
    Returns a list of directories to be used for downloading the files
    
    Args:
        years (list): The list of years
        months (list): The list of months
        days (list): The list of days
        level(list): The list of levels
        data_path (str): The path to the data folder (default: '../data/')
            This is the path relative to the current working directory
    
    Returns:
        dirs (list): The list of directories to be downloaded
    """

    # Create directories for the files
    # change to data folder from cwd
    create_directory(data_path)
    os.chdir(data_path)
    path = os.getcwd()
        
    # Create list of directories
    dirs = []
    for year in years:
        for month in months:
            for day in days:
                for _ in fileType:
                    for level in levels:
                        # dirs.append(path + '/Y' + year + '/M' + month + '/D' + day + '/')
                        dirs.append(path + '\\' + level + '\\Y' + year + '\\M' + month + '\\D' + day + '\\')

    return dirs

def buildURLList(years, months, days, fileTypes, levels):
    """
    Returns a list of urls to be downloaded based upon the user input
    
    Args:
        years (list): The list of years
        months (list): The list of months
        days (list): The list of days
        fileTypes (list): The list of file types
        levels (list): The list of levels
    
    Returns:
        urls (list): The list of urls to be downloaded
    """

    # MPLNET website
    main = 'https://mplnet.gsfc.nasa.gov/out/data/V3_partners/Appalachian_State/'

    # Create list of urls
    urls = []
    for year in years:
        for month in months:
            for day in days:
                for fileType in fileTypes:
                    for level in levels:
                        urls.append(main + 'Y' + year + '/M' + month + '/D' + day + '/MPLNET_V3_' + level + '_' + fileType + '_' + year + month + day + '_MPL44201' + '_Appalachian_State.nc4')

    return urls

def export(filename, files, variable):
    """
    Export data to csv file
    Uses current data directory to save csv file
    
    Args:
        filename (str): name of output csv file
        files (list): list of full file paths
        variable (str): name of variable to export
    
    Returns:
        None
    """

    # test if all files in files exist
    # if not, remove from list
    files = [x for x in files if os.path.isfile(x)]

    # handle case if file already exists
    if os.path.isfile(filename):
        # remove file
        os.remove(filename)

    # create empty dataframe
    df = pd.DataFrame()

    # loop through files
    # TODO: Error handling for variables that do not include time
    for file in files:
        try:
            # open file
            f = nc.Dataset(file, 'r')
            # get time
            time = f.variables['time'][:]
            # get altitude
            altitude = f.variables['altitude'][:]

            # get data
            data = None
            # Omit wavelength data as it is static, always 532 nm (green)
            # Check if wavelength data is present
            dims = list(f.variables[variable].dimensions)
            if len(dims) == 3:
                for num, dim in enumerate(dims):
                    if dim == 'wavelength':
                        if num == 0:
                            data = f.variables[variable][0, :, :]
                            dims = dims[1:]
                        elif num == 1:
                            data = f.variables[variable][:, 0, :]
                            dims = dims[0] + dims[2:]
                        elif num == 2:
                            data = f.variables[variable][:, :, 0]
                            dims = dims[:-1]
                    else:
                        if data is None:
                            data = f.variables[variable][:]
            elif len(dims) == 2:
                for num, dim in enumerate(dims):
                    if dim == 'wavelength':
                        if num == 0:
                            data = f.variables[variable][0, :]
                            dims.pop(0)
                        if num == 1:
                            data = f.variables[variable][:, 0]
                            dims.pop(1)
                    else: 
                        if data is None:
                            data = f.variables[variable][:]
            else: # only 1 diminsion
                data = f.variables[variable][:]

            # column naming
            if (len(data.shape) == 2):
                column_names = [str(x) for x in range(data.shape[1])]
            else:
                column_names = [variable]

            # add to dataframe
            if (df.empty == True):
                # only assign time if time if time is present
                if 'time' == dims[0]:
                    df = pd.DataFrame(data=data, index=time, columns=column_names)
                else:
                    df = pd.DataFrame(data=data, columns=column_names)
            else:
                if 'time' == dims[0]:
                    temp = pd.DataFrame(data=data, index=time, columns=column_names)
                    df = pd.concat([df, temp])
                else:
                    temp = pd.DataFrame(data=data, columns=column_names)
                    df = pd.concat([df, temp])

        except Exception as e:
            print(e)
            print('Error exporting file: {}'.format(file))
            continue
        finally:
            f.close()

    # rename index
    if dims[0] == 'time':
        df.index.rename('time', inplace=True)
    # if dims[1] is time flip data
    # TODO more testing
    if len(dims) > 1 and dims[1] == 'time':
        df = df.transpose()
        df.set_index(time, inplace=True)
        df.index.rename('time', inplace=True)


    # export to csv
    df.to_csv(filename, header=True)

def create_export_name(selectedVars, variable):
    """
    Create name for export file
    
    Args:
        selectedVars (object): object containing selected variables
        variable (str): name of variable to export

    Returns:
        filename (str): name of export file
    
    TODO: Include file type in filename
    """

    # get year range
    if len(selectedVars.selectedYears) > 1:
        year_range = 'Y{}-{}'.format(selectedVars.selectedYears[0], selectedVars.selectedYears[-1])
    else:
        year_range = 'Y' + selectedVars.selectedYears[0]

    # get month range
    if len(selectedVars.selectedMonths) > 1:
        month_range = 'M{}-{}'.format(selectedVars.selectedMonths[0], selectedVars.selectedMonths[-1])
    else:
        month_range = 'M' + selectedVars.selectedMonths[0]

    # get day range
    if len(selectedVars.selectedDays) > 1:
        day_range = 'D{}-{}'.format(selectedVars.selectedDays[0], selectedVars.selectedDays[-1])
    else:
        day_range = 'D' + selectedVars.selectedDays[0]

    level_range = selectedVars.selectedLevel[0]
    # create filename
    filename = '{}_{}_{}_{}_{}.csv'.format(variable.upper(), level_range, year_range, month_range, day_range)

    return filename

def returnDF(filename, files, variable):
    """
    Returns pandas DataFrame wtih data from all files
    
    Args:
        filename (str): name of output csv file
        files (list): list of full file paths
        variable (str): name of variable to export
    
    Returns:
        Pandas DataFrame
    """

    # test if all files in files exist
    # if not, remove from list
    files = [x for x in files if os.path.isfile(x)]

    # handle case if file already exists
    if os.path.isfile(filename):
        # remove file
        os.remove(filename)

    # create empty dataframe
    df = pd.DataFrame()

    # loop through files
    # TODO: Error handling for variables that do not include time
    for file in files:
        try:
            # open file
            f = nc.Dataset(file, 'r')
            # get time
            time = f.variables['time'][:]
            # get altitude
            altitude = f.variables['altitude'][:]

            # get data
            data = None
            # Omit wavelength data as it is static, always 532 nm (green)
            # Check if wavelength data is present
            dims = list(f.variables[variable].dimensions)
            if len(dims) == 3:
                for num, dim in enumerate(dims):
                    if dim == 'wavelength':
                        if num == 0:
                            data = f.variables[variable][0, :, :]
                            dims = dims[1:]
                        elif num == 1:
                            data = f.variables[variable][:, 0, :]
                            dims = dims[0] + dims[2:]
                        elif num == 2:
                            data = f.variables[variable][:, :, 0]
                            dims = dims[:-1]
                    else:
                        if data is None:
                            data = f.variables[variable][:]
            elif len(dims) == 2:
                for num, dim in enumerate(dims):
                    if dim == 'wavelength':
                        if num == 0:
                            data = f.variables[variable][0, :]
                            dims.pop(0)
                        if num == 1:
                            data = f.variables[variable][:, 0]
                            dims.pop(1)
                    else: 
                        if data is None:
                            data = f.variables[variable][:]
            else: # only 1 diminsion
                data = f.variables[variable][:]

            # column naming
            if (len(data.shape) == 2):
                column_names = [str(x) for x in range(data.shape[1])]
            else:
                column_names = [variable]

            # add to dataframe
            if (df.empty == True):
                # only assign time if time if time is present
                if 'time' == dims[0]:
                    df = pd.DataFrame(data=data, index=time, columns=column_names)
                else:
                    df = pd.DataFrame(data=data, columns=column_names)
            else:
                if 'time' == dims[0]:
                    temp = pd.DataFrame(data=data, index=time, columns=column_names)
                    df = pd.concat([df, temp])
                else:
                    temp = pd.DataFrame(data=data, columns=column_names)
                    df = pd.concat([df, temp])

        except Exception as e:
            print(e)
            print('Error exporting file: {}'.format(file))
            continue
        finally:
            f.close()

    # rename index
    if dims[0] == 'time':
        df.index.rename('time', inplace=True)
    # if dims[1] is time flip data
    # TODO more testing
    if len(dims) > 1 and dims[1] == 'time':
        df = df.transpose()
        df.set_index(time, inplace=True)
        df.index.rename('time', inplace=True)

    # return the dataframe
    return df

def MinuteAvg(filename, files, variable):
    """
    Returns initial dataframe formatted the same as the hr and day avg dataframes as data is already in minute bins.
    If the data is 3D/altitude based the column headers are the altitudes in km.
    If the data is 2D then the column header is the variable.  
    
    Args:
        filename (str): name of output csv file of the data
        files (str): list of full file paths
        variable (str): name of the variable selected

    Returns:
        Pandas DataFrame
    """
    #create a dataframe of the data
    df = returnDF(filename, files, variable)

    #Makes an array of the altitudes 
    alt_arr = np.linspace(1.1549481, 31.059248, 400)
    alt_arr_str = ['']*400

    #create array of Julian Day Numbers incremented by 1 hour starting at Midnight of the first selected day ending at 23:00 of final day
    time_min = np.linspace(df.index[0]-(30/(24*60*60)), df.index[len(df.index)-1]-(30/(24*60*60)), int((len(df.index))))
    time_min = pd.to_datetime(time_min, unit = 'D', origin = 'julian').round(freq='T')
    if (variable == 'cloud_base' or variable == 'cloud_top'): 
        #if cld_base or cld_top, change to a 1 column array of 1s and 0s, 1 if there is a cloud and 0 if there is not a cloud.
        df = df.drop(df.iloc[:, 1:], axis=1)
        df.loc[~df['0'].isna(), '0'] = 1
        df.loc[df['0'].isna(), '0'] = 0
        df.columns = [variable]
    elif len(df.columns) == 400:
        #if the dataframe is instead 400 columns this indicates 3D data with altitude bins at each time.
        #This just makes an array of altitude strings for the column headers
        for i in range(len(alt_arr)):
            alt_arr_str[i] = '{:.3f}'.format(alt_arr[i])+' km'
        df.columns = alt_arr_str
        
    df.set_index(time_min, inplace=True) #sets the index to the time
    df.index.rename('time', inplace=True) #renames index header to "time"
    df.dropna(how='all', axis = 'columns', inplace = True) #drops all empty columns
    return df
    
def HrAvg(filename, files, variable):
    """
    Returns an hourly averaged pandas DataFrame with the time as the index.
    If the data is 3D/altitude based the column headers are the altitudes in km.
    If the data is 2D then the column header is the variable. 
    
    Args:
        filename (str): name of output csv file of the data
        files (str): list of full file paths
        variable (str): name of the variable selected

    Returns:
        Pandas DataFrame
    """
    #create a dataframe of the data
    df = returnDF(filename, files, variable)

    #Makes an array of the altitudes
    alt_arr = np.linspace(1.1549481, 31.059248, 400)
    alt_arr_str = ['']*400

    #create array of Julian Day Numbers incremented by 1 hour starting at Midnight of the first selected day ending at 23:00 of final day
    #might not need this line as the resample function might do it for you
    time_hrs = np.linspace(df.index[0]-(30/(24*60*60)), df.index[len(df.index)-60]-(30/(24*60*60)), int((len(df.index))))
    df['Datetime'] = pd.to_datetime(time_hrs, unit = 'D', origin = 'julian').round(freq='h')

    #makes the index of the dataframe the datetimes
    df = df.set_index('Datetime')
    if (variable == 'cloud_base' or variable == 'cloud_top'):
        #if cld_base or cld_top, change to a 1 column array of 1s and 0s, 1 if there is a cloud and 0 if there is not a cloud. 
        #hourly average will then show cloud density for that hour from 0-1
        df = df.drop(df.iloc[:, 1:], axis=1)
        df.loc[~df['0'].isna(), '0'] = 1
        df.loc[df['0'].isna(), '0'] = 0
        df.columns = [variable] #makes column header cloud base or cloud top
    elif len(df.columns) == 400:
        #if the dataframe is instead 400 columns this indicates 3D data with altitude bins at each time.
        #This just makes an array of altitude strings for the column headers
        for i in range(len(alt_arr)):
            alt_arr_str[i] = '{:.3f}'.format(alt_arr[i])+' km'
        df.columns = alt_arr_str
    else:
        #if neither of these things indicating non altitude based data and therefore 2d data, make the column header the variable
        df.columns = [variable]

    df = df.resample('h').mean() #take an hourly mean of the whole dataframe

    df.index.rename('time', inplace=True) #renames the df index header to time
    df.dropna(how='all', axis = 'columns', inplace = True) #drops all empty columns
    return df

def DayAvg(filename, files, variable):
    """
    Returns a daily averaged pandas DataFramewith the time as the index in pandas DateTime format.
    If the data is 3D/altitude based the column headers are the altitudes in km.
    If the data is 2D then the column header is the variable.
    
    Args:
        filename (str): name of output csv file of the data
        files (str): list of full file paths
        variable (str): name of the variable selected

    Returns:
        Pandas DataFrame
    """
    #create a dataframe of the data
    df = returnDF(filename, files, variable)

    #create an array of the altitudes to use as the headers of the final dataframe
    alt_arr = np.linspace(1.1549481, 31.059248, 400)
    alt_arr_str = ['']*400

    #converts the julian day times to pandas datetime objects
    df['Datetime'] = pd.to_datetime(df.index, unit = 'D', origin = 'julian').round(freq='s')

    #create a dataframe with the first column of the dataframe as the datetimes
    df = df.set_index('Datetime')
    if (variable == 'cloud_base' or variable == 'cloud_top'):
        #if cld_base or cld_top, change to a 1 column array of 1s and 0s, 1 if there is a cloud and 0 if there is not a cloud. 
        #daily average will then show cloud density for that day from 0-1
        df = df.drop(df.iloc[:, 1:], axis=1)
        df.loc[~df['0'].isna(), '0'] = 1
        df.loc[df['0'].isna(), '0'] = 0
        df.columns = [variable]
    elif len(df.columns) == 400:
        #if the dataframe is instead 400 columns this indicates 3D data with altitude bins at each time.
        #This just makes an array of altitude strings for the column headers
        for i in range(len(alt_arr)):
            alt_arr_str[i] = '{:.3f}'.format(alt_arr[i])+' km'
        df.columns = alt_arr_str
    else:
        #if neither of these things indicating non altitude based data and therefore 2d data, make the column header the variable
        df.columns = [variable]

    df = df.resample('D').mean() #take a daily average of the whole dataframe
     
    df.index.rename('time', inplace=True) #rename the index header to 'time'
    df.dropna(how='all', axis = 'columns', inplace = True) #drop all empty columns
    return df