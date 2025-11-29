from ebas.io.file import nasa_ames
from nilutility.datatypes import DataObject
from ebas.domain.basic_domain_logic.time_period import estimate_period_code, \
    estimate_resolution_code, estimate_sample_duration_code
import datetime
import pandas as pd 
import numpy as np 
from pathlib import Path

__version__ = '1.00.00'

def set_fileglobal_metadata(nas,name,email,instrument):
    """
    Set file global metadata for the EbasNasaAmes file object

    Parameters:
        nas: EbasNasaAmes file object
        name: name as a list of str ['first', 'last']
        email: email as a str 
        instrument: intsrument and actris type as a list of str ['instrument', 'GAW instrument type']
    Returns:
        None
    """
    # All times reported to EBAS need to be in UTC!
    # Setting the timezone here explicitly should remind you to check your data
    nas.metadata.timezone = 'UTC'

    # Revision information
    nas.metadata.revdate = datetime.datetime.now()
    nas.metadata.revision = '1.1a'
    nas.metadata.revdesc = \
        'initiol revision to ebas, generated with MyDataTool 1.22'

    # Data Originator Organisation
    nas.metadata.org = DataObject(
        OR_CODE='APP', #<--- ask jps
        OR_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
        OR_ACRONYM='AppalAIR', OR_UNIT='Department of Physics & Astronomy',
        OR_ADDR_LINE1='525 Rivers Street', OR_ADDR_LINE2=None,
        OR_ADDR_ZIP='28608', OR_ADDR_CITY='Boone', OR_ADDR_COUNTRY='United States of America')

    # Projects the data are associated to
    nas.metadata.projects = ['GAW-WDCA', 'AppalAIR']

    # Data Originators (PIs)
    nas.metadata.originator = []
    nas.metadata.originator.append(
        DataObject(
            PS_LAST_NAME='James', PS_FIRST_NAME='Sherman', PS_EMAIL='shermanjp@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics & Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))
    nas.metadata.originator.append(
        DataObject(
            PS_LAST_NAME=name[-1], PS_FIRST_NAME=name[0], PS_EMAIL=email,
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics & Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))

    # Data Submitters (contact for data technical issues)
    nas.metadata.submitter = []
    nas.metadata.submitter.append(
        DataObject(
            PS_LAST_NAME=name[-1], PS_FIRST_NAME=name[0], PS_EMAIL=email,
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics & Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))

    # Station metadata
    nas.metadata.station_code = 'APP' #<--ASK Dr. Sherman
    nas.metadata.platform_code = 'APP'
    nas.metadata.station_name = u'AppalAIR'

    nas.metadata.station_wdca_id = '2'#<--ASK Dr. Sherman'
    nas.metadata.station_gaw_id = '2'#<--ASK Dr. Sherman
    nas.metadata.station_gaw_name = u'AppalAIR'
    # nas.metadata.station_airs_id =    # N/A
    nas.metadata.station_other_ids = '2'#<--ASK Dr. Sherman'
    # nas.metadata.station_state_code =  # N/A
    nas.metadata.station_landuse = 'Residential'
    nas.metadata.station_setting = 'Mountain'
    nas.metadata.station_gaw_type = '2'#<--ASK Dr. Sherman'
    nas.metadata.station_wmo_region = 4
    nas.metadata.station_latitude = 36.212801
    nas.metadata.station_longitude = -81.692592
    nas.metadata.station_altitude = 1079

    # More file global metadata, but those can be overridden per variable
    # See set_variables for examples
    nas.metadata.instr_type = instrument[-1]
    nas.metadata.lab_code = '2'#<--ASK Dr. Sherman'
    nas.metadata.instr_name = instrument[0]
    nas.metadata.method = f'LABCODE_{instrument[0]}'
    nas.metadata.regime = 'IMG'
    nas.metadata.matrix = 'aerosol'
    #nas.metadata.comp_name   will be set on variable level
    #nas.metadata.unit        will be set on variable level
    nas.metadata.statistics = 'arithmetic mean'
    nas.metadata.datalevel = '2'

def set_time_axes(nas, dates):
    """
    Set the time axes and related metadata for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
        dates  list of dates for datafile 
    Returns:
        None
    """
    # define start and end times for all samples
    date_list =[]
    for i in range(len(dates)):
        date_list.append([pd.to_datetime(dates[i]),pd.to_datetime(dates[i])+pd.Timedelta(59, 'min')])  
    # converts to a time stampt and appends the first minute and last minute of the sampled hour
    input(date_list)
    nas.sample_times = \
        date_list

    # Generate metadata that are related to the time axes:

    # period code is an estimate of the current submissions period, so it should
    # always be calculated from the actual time axes, like this:
    nas.metadata.period = estimate_period_code(nas.sample_times[0][0],
                                               nas.sample_times[-1][1])

    # Sample duration can be set automatically
    nas.metadata.duration = estimate_sample_duration_code(nas.sample_times)
    # or set it hardcoded:
    # nas.metadata.duration = '3mo'

    # Resolution code can be set automatically
    # But be aware that resolution code is an identifying metadata element.
    # That means, several submissions of data (multiple years) will
    # only be stored as the same dataset if the resolution code is the same
    # for all submissions!
    # That might be a problem for time series with varying resolution code
    # (sometimes 2 months, sometimes 3 months, sometimes 9 weeks, ...). You
    # might consider using a fixed resolution code for those time series.
    # Automatic calculation (will work from ebas.io V.3.0.7):
    nas.metadata.resolution = estimate_resolution_code(nas.sample_times)
    # or set it hardcoded:
    # nas.metadata.resolution = '3mo'

    # It's a good practice to use Jan 1st of the year of the first sample
    # endtime as the file reference date (zero point of time axes).
    nas.metadata.reference_date = \
        datetime.datetime(nas.sample_times[0][1].year, 1, 1)

def set_variables(nas,data_table, flag_table, headers):
    """
    Set metadata and data for all variables for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
        data_table:  list of lists for data. Each row represents a csv column
        flag_table  list of list of lists. Each row of list represents a csv column where each 
                    list within the row is the total flags for the ebas value since each row can 
                    have multiple values
        headers     list of columnn names for out putted date file
    Returns:
        None
    """
    for i in range(len(data_table)): # flags and headers need to be same shape as data
        values = data_table[i]
        flags = flag_table[i]
        metadata = DataObject()
        metadata.comp_name = headers[i]
        metadata.unit = '1/cm3'
        nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,
                                    metadata=metadata))

def ebas_genfile(path,data,flag,date, headers, name, email, inst):
    """
    Main program for ebas_flatcsv
    Created for lexical scoping.

    Inputs:
        path: Folder path for outputs
        data: Data for processing. List of lists such that rows represents original
          columns of data
        flag: Flags to correspond to data. List of list of lists. Each row is a column to 
          append to data, and each row is made of of lists of 'flags' such that each each 
          row can have multiple flags
        date: Dates for data. 
        headers: Column names for each row of data
        name: list for name of orriginator ['first','last']
        email: email for originator str
        inst: list for instrument data ['instrument name', 'GAW instrument type']

    Returns:
        none
    """

    # Create an EbasNasaAmes file object
    nas = nasa_ames.EbasNasaAmes()

    # Set file global metadata
    set_fileglobal_metadata(nas,name,email, inst)

    # Set the time axes and related metadata
    set_time_axes(nas,date)
    # Set metadata and data for all variables
    set_variables(nas, data, flag, headers)

    # write the file:
    nas.write(createfiles=True, destdir=path)
    # createfiles=True
    #     Actually creates output files, else the output would go to STDOUT.
    # You can also specify:
    #     destdir='path/to/directory'
    #         Specify a specific relative or absolute path to a directory the
    #         files should be written to
    #     flags=FLAGS_COMPRESS
    #         Compresses the file size by reducing flag columns.
    #         Flag columns will be less explicit and thus less intuitive for
    #         humans to read.
    #     flags=FLAGS_ALL
    #         Always generate one flag column per variable. Very intuitive to
    #         read, but increases filesize.
    #     The default for flags is: Generate one flag column per file if the
    #     flags are the same for all variables in the file. Else generate one
    #     flag column per variable.
    #     This is a trade-off between the advantages and disadvantages of the
    #     above mentioned approaches.

def quickfill(csv):
    '''Self contained method for formatting csv data into NASA ames files
    without having to know any python. Simply pass a formated CSV and respond
    to inputs in the terminal:
    Inputs:
        csv: path to comma seperated values file to process. 
            FORMAT REQUIREMENTS:
            - hourly averaged 
            - first column is datetime index
            - ACTRIS flag column in csv
            - columns are identical between csv
              and final actris formating
    Returns:
        none
    '''
    df = pd.read_csv(csv)#read in csv 
    df.index = pd.to_datetime(df.index) 
    dates= df.index.to_list()
    cols = input('list of column names seperated by a comma: ').split(',')
    data = df[cols].values.tolist()
    if input('is there a flag column for the data?(Y/N)') == 'Y':
        flag_col = input('flag column header: ')
        flags = df[flag_col].values.tolist()
    else:
        flags = [[[000] for j in range(len(data[i]))] for i in range(len(dates))]
    data =  list(map(list, zip(*data))) #Transpose data
    flags = list(map(list, zip(*flags))) #Transpose flags
    name = input('full name ["first last"]:').split()
    email = input('email address')
    inst = input('instrument name and GAW instrument type seperated by a comma:').split(',')
    path = Path(input('path to output directory: '))
    ebas_genfile(path,data,flags,dates, cols, name, email, inst)

quickfill(Path(input('Path to csv:')))