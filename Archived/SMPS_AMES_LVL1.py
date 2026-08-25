#!/usr/bin/env python
# coding=utf-8
"""
$Id: ebas_genfile.py 1431 2016-10-27 14:25:12Z pe $

Example for creating an EBAS_1.1 NasaAmes datafile.
"""
from ebas.io.file import nasa_ames
from nilutility.datatypes import DataObject
from ebas.domain.basic_domain_logic.time_period import estimate_period_code, \
    estimate_resolution_code, estimate_sample_duration_code
import datetime
from ebas.io.ebasmetadata import DatasetCharacteristicList
import pandas as pd
import numpy
def testnan(val):
    if val * 0. == 0:
        return val
    else:
        return None

#INPUT DATA-VALUES
f = r'Level1_final.txt'
first_column=1 #(first column with variable to process) 
last_column=75 # (last column to process)
index_flag_col=76 #(nb column with flag data)
tresol=5 #sample duration
df=pd.read_csv(f, delimiter=',', header=0, engine='python')  #Specify nb_lines of header, delimiter etc
df = df.where((pd.notnull(df)), None) #Replace nan values with None
start_times_df=pd.to_datetime(df['DateTime'],format='%Y-%m-%d %H:%M:%S') #Specify format of datetime column to read. The column in called by the header name
names=list(df.columns.values) 
numb_col=df.shape[1]

#INPUT METADATA_VARIABLES !!!! METADATA and INPUT FILES MUST HAVE THE SAME NUMBER OF COLUMNS AND IN THE SAME ORDER AND POSITION!!
f_metadata=r'Input_metadata_lev1.csv'
df_metadata=pd.read_csv(f_metadata, delimiter=',', header=0)

#OUTPUT FOLDER
destdir_out='./'




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
    nas.metadata.timezone = 'UTC'
    
    nas.metadata.revdate = datetime.datetime(2017, 7, 20, 13, 23, 00)
    nas.metadata.revision = '1.1a'
    nas.metadata.revdesc = \
        'initial revision, NRT SMPS_control_lev0 lev0 2 lev1 v.0.0_2'
   
    nas.metadata.type = 'TU'
  
    # Revision information
    nas.metadata.startdate = datetime.datetime(2016, 1, 1, 00, 00, 00)#not ready #I am not sure about the time format


    nas.metadata.datalevel = '1'
    nas.metadata.period = '1y'
    nas.metadata.resolution = '5mn'
    nas.metadata.duration = '5mn'
    nas.metadata.rescode_sample = '5mn'
    
    # Data Originator Organisation
    nas.metadata.org = DataObject(
            OR_CODE='APP', #not ready#NOT RIGHT, ASK JPS
            OR_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            OR_ACRONYM='AppalAIR', OR_UNIT='Department of Physics and Astronomy',
            OR_ADDR_LINE1='525 Rivers Street', OR_ADDR_LINE2=None,
            OR_ADDR_ZIP='28608', OR_ADDR_CITY='Boone', OR_ADDR_COUNTRY='United States of America')


    # Projects the data are associated to
    nas.metadata.projects = ['GAW-WDCA', 'AppalAIR']
# Station metadata
    nas.metadata.station_code = 'NO0042G' #not ready#ASK JPS FOR STATION AND PLATFORM CODES
    nas.metadata.platform_code = 'NO0042S'#not ready
    nas.metadata.station_name = u'AppalAIR'

    nas.metadata.station_wdca_id = 'GAWANO__ZEP'#not ready
    nas.metadata.station_gaw_id = 'ZEP'#not ready
    nas.metadata.station_gaw_name = u'AppalAIR'
    # nas.metadata.station_airs_id =    # N/A
    nas.metadata.station_other_ids = '721 (NILUDB)'#not ready
    # nas.metadata.station_state_code =  # N/A
    nas.metadata.station_landuse = 'Residential'
    nas.metadata.station_setting = 'Mountain'
    nas.metadata.station_gaw_type = 'G'
    nas.metadata.station_wmo_region = 4
    nas.metadata.station_latitude = 36.212801
    nas.metadata.station_longitude = -81.692592
    nas.metadata.station_altitude = 1079.0
    
    # More file global metadata, but those can be overridden per variable
    # See set_variables for examples
    
    nas.metadata.comp_name = 'particle_number_size_distribution'
    nas.metadata.unit = '1/cm3'#not ready
    nas.metadata.matrix = ''#not ready
    nas.metadata.lab_code = ''#not ready
    nas.metadata.instr_type = 'smps'#not ready #MUST BE in LOWER CASE LETTERS!! OTHERWISE GIVES ERRORS
    nas.metadata.instr_name = ''#not ready
    nas.metadata.instr_manufacturer ='TROPOS-TSI'#not ready
    nas.metadata.instr_model = 'TROPOS-SMPS'#not ready
    nas.metadata.instr_serialno = 'SMPS=TROPOS' #not ready
    
    nas.metadata.method = 'GR05L__NRT_SMPS_lev1'#not ready
    nas.std_method = ''#not ready
    
    nas.metadata.inlet_type = 'Impactor direct'#not ready
    nas.metadata.inlet_desc = 'PM10 impactor'#not ready
    nas.metadata.hum_temp_ctrl = 'Nafion dryer'#not ready
    nas.metadata.hum_temp_ctrl_desc = 'Humidity/temperature control description: sample dried to below 40% RH with nafion dryer'#not ready
    
    nas.metadata.vol_std_temp = 273.15 #not ready
    nas.metadata.vol_std_pressure = 1013.25  #not ready
    nas.metadata.detection_limit=(0, '1/cm3') #not ready
    nas.metadata.detection_limit_desc='Determined only by instrument counting statistics, no detection limit flag used'#not ready
    nas.metadata.uncertainty_desc='uncertainty range between instruments in intercomparison by Jiang et al. 2014.'#not ready
   
    nas.metadata.zero_negative='Zero values possible'#not ready
    nas.metadata.zero_negative_desc='Zero values may appear due to statistical variations at very low concentrations'#not ready
   
    
 
    # Data Originators (PIs)
    nas.metadata.originator = []
    nas.metadata.originator.append(
        DataObject(
            PS_LAST_NAME='Sherman', PS_FIRST_NAME='James', PS_EMAIL='shermanjp@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics and Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))
    nas.metadata.originator.append(
        DataObject(
            PS_LAST_NAME=u'Parkhurst', PS_FIRST_NAME='Ethan', PS_EMAIL='parkhursted@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics and Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))
    nas.metadata.originator.append(
        DataObject(
            PS_LAST_NAME=u'Kitteringham', PS_FIRST_NAME='Ryan', PS_EMAIL='kitteringhamrr@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics and Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))

    # Data Submitters (contact for data technical issues)
    nas.metadata.submitter = []
    nas.metadata.submitter.append(
        DataObject(
            PS_LAST_NAME=u'Parkhurst', PS_FIRST_NAME='Ethan', PS_EMAIL='parkhursted@appstate.edu',
            PS_ORG_NAME='Appalachian Atmospheric Interdisciplinary Research Program',
            PS_ORG_ACR='AppalAIR', PS_ORG_UNIT='Department of Physics and Astronomy',
            PS_ADDR_LINE1='525 Rivers Street', PS_ADDR_LINE2=None,
            PS_ADDR_ZIP='28608', PS_ADDR_CITY='Boone',
            PS_ADDR_COUNTRY='United States of America',
            PS_ORCID=None,
        ))

       
    #nas.metadata.qm_id=''
    #nas.metadata.qm_doc_date=''
    #nas.metadata.qm_doc_url=''
    
    nas.metadata.acknowledgements='Request acknowledgment details from data originator'#not ready
    nas.metadata.comment='none'#not ready

def set_time_axes(nas):
    """
    Set the time axes and related metadata for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
    Returns:
        None
    """
    # define start and end times for all samples
    start_times_list=start_times_df.tolist()
    start_times_dt_list=[start_times_list[i].to_pydatetime()for i in range(len(start_times_list))]
    end_times_df=start_times_df+datetime.timedelta(minutes=tresol)
    end_times_list=end_times_df.tolist()
    end_times_dt_list=[end_times_list[i].to_pydatetime()for i in range(len(start_times_list))]
    nas.sample_times=list(zip(start_times_dt_list,end_times_dt_list))
    #nas.sample_times = \
    #    [(datetime.datetime(2008, 1, 1, 0, 0), datetime.datetime(2008, 1, 1, 1, 0)),
    #     (datetime.datetime(2008, 1, 1, 1, 0), datetime.datetime(2008, 1, 1, 2, 0)),
    #     (datetime.datetime(2008, 1, 1, 2, 0), datetime.datetime(2008, 1, 1, 3, 0)),
    #     (datetime.datetime(2008, 1, 1, 3, 0), datetime.datetime(2008, 1, 1, 4, 0)),
    #     (datetime.datetime(2008, 1, 1, 4, 0), datetime.datetime(2008, 1, 1, 5, 0)),
    #     (datetime.datetime(2008, 1, 1, 5, 0), datetime.datetime(2008, 1, 1, 6, 0)),
    #     (datetime.datetime(2008, 1, 1, 6, 0), datetime.datetime(2008, 1, 1, 7, 0)),
    #     (datetime.datetime(2008, 1, 1, 7, 0), datetime.datetime(2008, 1, 1, 8, 0))]

    #
    # Generate metadata that are related to the time axes:
    #

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

def set_variables(nas):
    """
    Set metadata and data for all variables for the EbasNasaAmes file object.

    Parameters:
        nas    EbasNasaAmes file object
    Returns:
        None
    """
    
    #x=14
    x=first_column
    while x < last_column+1:                                                  
        col=df.take([x],axis=1)
        values=col.values.T.tolist()[0]
        values=numpy.array(values, dtype=numpy.float)
        values=[ round(elem,6) for elem in values]
        values=list(map(testnan, values))
        col_metadata=df_metadata.take([x],axis=1)
        col_flag=df.take([index_flag_col],axis=1).astype(int)
        flags=col_flag.values.tolist() 
        metadata = DataObject()
        metadata.comp_name = col_metadata[0:1].values.T.tolist()[0][0]
       
        #Unit metadata
        metadata.unit = col_metadata[1:2].values.T.tolist()[0][0]
        if metadata.unit!=metadata.unit:
            metadata.unit = None
        else:
            metadata.unit = metadata.unit
                  
        metadata.matrix=col_metadata[2:3].values.T.tolist()[0][0]
        metadata.title=col_metadata[3:4].values.T.tolist()[0][0]   
        
        #Statistics metadata
        metadata.statistics=col_metadata[4:5].values.T.tolist()[0][0]
        if metadata.statistics!=metadata.statistics:
            metadata.statistics = None
        else:
            metadata.statistics = metadata.statistics
            
        #Uncertainity metadata    
        uncertainty_value=col_metadata[5:6].values.T.tolist()[0][0]
        uncertainty_unit=col_metadata[6:7].values.T.tolist()[0][0]
        if uncertainty_value!=uncertainty_value:
            metadata.uncertainty = None
        else:
            metadata.uncertainty = (int(uncertainty_value), uncertainty_unit)
            
        

        #Characteristics    
        if col_metadata[7:8].values.T.tolist()[0][0]=='None':
            nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,
                                    metadata=metadata))
        else:
            metadata.characteristics = DatasetCharacteristicList()
            Characteristic_type=col_metadata[7:8].values.T.tolist()[0][0]
            Value=col_metadata[8:9].values.T.tolist()[0][0]
            Instrument_type=col_metadata[9:10].values.T.tolist()[0][0]
            metadata.characteristics.add_parse(Characteristic_type, Value, Instrument_type, metadata.comp_name)  
            nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,
                                metadata=metadata))

            
        #nas.variables.append(DataObject(values_=values, flags=flags, flagcol=True,
                                    #metadata=metadata))
                                    #
        x+=1    

def ebas_genfile():
    """
    Main program for ebas_flatcsv
    Created for lexical scoping.

    Parameters:
        None
    Returns:
        none
    """

    # Create an EbasNasaAmes file object
    nas = nasa_ames.EbasNasaAmes()

    # Set file global metadata
    set_fileglobal_metadata(nas)

    # Set the time axes and related metadata
    set_time_axes(nas)

    # Set metadata and data for all variables
    set_variables(nas)

    # write the file:
    nas.write(createfiles=True,destdir=destdir_out)
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

ebas_genfile()
