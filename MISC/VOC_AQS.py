import os
import time
import json
from os.path import expanduser 

import requests
import pandas as pd
from download_aqs import main as pullAQI



OUTDIR_JSON = expanduser("~/Documents/Research/AQS_Json/")