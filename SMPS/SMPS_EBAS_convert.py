from ebas.io.file import nasa_ames
from nilutility.datatypes import DataObject
from ebas.domain.basic_domain_logic.time_period import estimate_period_code, \
    estimate_resolution_code, estimate_sample_duration_code
import datetime
import pandas as pd 
import numpy as np 
from pathlib import Path
import os

__version__ = '1.00.00'

def set_fileglobal_metadata(nas):
    """
    Set file global metadata for the EbasNasaAmes file object

    Parameters:
        nas    EbasNasaAmes file object
    Returns:
        None
    """
    # All times reported to EBAS need to be in UTC!
    # Setting the timezone here explicitly should remind you to check your data
    nas.metadata.typecode = 'TI'   # time irregular-> time series with gaps
    nas.metadata.timezone = 'UTC'

    # Revision information
    nas.metadata.revdate = datetime.datetime.now()
    nas.metadata.revision = '1.1a'
    nas.metadata.revdesc = \
        'initial revision to ebas, generated with MyDataTool 1.22'

    # Data Originator Organisation
    nas.metadata.org = DataObject(
        OR_CODE='US10L',
        OR_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
        OR_ACRONYM='AppalAIR', OR_UNIT='Department of Physics & Astronomy',
        OR_ADDR_LINE1='525 Rivers Street', OR_ADDR_LINE2=None,
        OR_ADDR_ZIP='28608', OR_ADDR_CITY='Boone', OR_ADDR_COUNTRY='United States of America')

    # Projects the data are associated to
    nas.metadata.projects = ['GAW-WDCA', 'NOAA-ESRL']

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
            PS_LAST_NAME=u'Sykes', PS_FIRST_NAME='Ben', PS_EMAIL='sykesbb@appstate.edu',
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
            PS_LAST_NAME=u'Sykes', PS_FIRST_NAME='Ben', PS_EMAIL='sykesbb@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics & Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))

    # Station metadata
    nas.metadata.station_code = 'US3446S' 
    nas.metadata.platform_code = 'US3446C'
    nas.metadata.station_name = u'AppalAIR'

    nas.metadata.station_gaw_id = 'APP'
    nas.metadata.station_gaw_name = u'Appalachian State University'
    # nas.metadata.station_airs_id =    # N/A
    # nas.metadata.station_other_ids = '2'
    nas.metadata.station_state_code =  'US3446C'
    nas.metadata.station_landuse = 'Residential'
    nas.metadata.station_setting = 'Mountain'
    nas.metadata.station_gaw_type = 'R'
    nas.metadata.station_wmo_region = 4
    nas.metadata.station_latitude = 36.2130012512
    nas.metadata.station_longitude = -81.6920013428
    nas.metadata.station_altitude = 1076 #[m]
    

    # More file global metadata, but those can be overridden per variable
    # See set_variables for examples
    nas.metadata.instr_type = 'smps'
    nas.metadata.lab_code = 'US10L'
    nas.metadata.instr_name = 'TSI_SMPS3938_APP'
    nas.metadata.instr_manufacturer ='TSI'
    nas.metadata.instr_model = '3938'
    nas.metadata.method = 'US10L_smps'
    nas.metadata.regime = 'IMG'
    nas.metadata.matrix = 'aerosol'
    nas.std_method = 'SOP=Wiedensohler2012' 
    
    nas.metadata.inlet_type = 'Open conductive tubing'
    nas.metadata.hum_temp_ctrl = 'Nafion dryer'
    nas.metadata.hum_temp_ctrl_desc = 'sample dried to below 40% RH with nafion dryeroperated with reduced pressure sheath air'
    
    # nas.metadata.comp_name ='particle_number_size_distribution'
    # nas.metadata.unit  ='#/cm3'    #will be set on variable level
    nas.metadata.vol_std_temp = 273.15 #[K]
    nas.metadata.vol_std_pressure = 1013.25 #[hpa]
    nas.metadata.statistics = 'arithmetic mean'
    nas.metadata.datalevel = '2'
    nas.metadata.period = '1y'
    nas.metadata.resolution = '5mn'
    nas.metadata.duration = '5mn'
    nas.metadata.rescode_sample = '5mn'
    nas.metadata.detection_limit=(0, '1/cm3') 
    nas.metadata.detection_limit_desc=' Determined only by instrument counting statistics, no detection limit flag used'
    nas.metadata.uncertainty_desc='uncertainty range between instruments in intercomparison by Wiedensohler et al. 2012. (AMT)'
   
    nas.metadata.zero_negative='Zero possible'
    nas.metadata.zero_negative_desc='Zero values may appear due to statistical variations at very low concentrations'
   
    

def set_time_axes(nas, dates):
    """
    Set the time axes and related metadata for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
    Returns:
        None
    """
    # define start and end times for all samples
    date_vals = []
    date_list =[]
    for i in range(len(dates)):
        date_list.append([pd.to_datetime(dates[i]),pd.to_datetime(dates[i])+pd.Timedelta(59, 'min')])
    #     if i == 0:
    #         date_vals.append(pd.to_datetime(dates[i]))
    #     elif (pd.to_datetime(dates[i+1])-pd.to_datetime(dates[i]) > pd.Timedelta(1, 'hr')):
    #         date_vals.append(pd.to_datetime(dates[i]))
    #         date_vals.append(pd.to_datetime(dates[i+1]))
    # date_vals.append(pd.to_datetime(dates[-1]))   
    # date_list = [date_vals[i:i + 2] for i in range(0, len(date_vals), 2)]           

    nas.sample_times = \
        date_list
        # [(datetime.datetime(2025, 1, 2, 15, 0), datetime.datetime(2025, 10, 31, 23, 0))]
    
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

# def set_variables(nas,data_table, flag_table, headers,cols):
#     """
#     Set metadata and data for all variables for the EbasNasaAmes file object.

#     Parameters:
#         nas    EbasNasaAmes file object
#     Returns:
#         None
#     """
#     # variable 1: examples for missing values and flagging
#     for i in range(len(data_table)):
#         values = data_table[i]
#         flags = flag_table[i]#[[0] for i in range(len(values))] #
#         # input(np.shape(values))
#         # values = np.nan_to_num(data, nan=0)
#         metadata = DataObject()
#         metadata.comp_name = headers[i]
#         # metadata.unit = 
#         metadata.title = cols[i]
#         nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,
#                                     metadata=metadata))

def set_variables(nas,data_table, flag_table, headers, units, stats, locs, Ds, cols):
    """
    Set metadata and data for all variables for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
    Returns:
        None
    """
    # variable 1: examples for missing values and flagging
    for i in range(len(data_table)):
        values = data_table[i]
        flags = flag_table[i]
        metadata = DataObject()
        metadata.comp_name = headers[i]
        metadata.unit = units[i]
        metadata.statistics = stats[i]
        metadata.title = cols[i]
        nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,metadata=metadata))
        if locs[i] != '':
            nas.add_var_characteristics(-1, 'Location', locs[i])
        if Ds[i] != '':
            nas.add_var_characteristics(-1, 'D', Ds[i])

def ebas_genfile(path, data, flag, date, headers, units, stats, locs, Ds, cols): #ebas_genfile(path,data,flag,date, headers,cols):
    """
    Main program for ebas_flatcsv
    Created for lexical scoping.

    Inputs:
        path: Folder path for outputs
        data: Data for processing. List of lists such that rows represents original
          columns of data
        flag: Flags to correspond to data. List of list of lists. Each row is a column to 
          append to data, and each row is made of of lists of 'flags' such that each 
    Returns:
        none
    """
    path_out = path

    # Create an EbasNasaAmes file object
    nas = nasa_ames.EbasNasaAmes()

    # Set file global metadata
    set_fileglobal_metadata(nas)

    # Set the time axes and related metadata
    set_time_axes(nas,date)
    # Set metadata and data for all variables
    set_variables(nas, data, flag, headers, units, stats, locs, Ds, cols)

    # write the file:
    # with open(r"C:\Users\bensy\Documents\Research\EBAStest.txt", "w") as text_file:
    #     text_file.write(str(nas.variables))
    nas.write(createfiles=True, destdir=path_out)
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
