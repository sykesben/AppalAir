### AppState MPLNet

MPLNet data tool for use with AppState MPLNet data.
This tool provides functionality to download MPLNet data, convert from NetCDF4 to CSV, and some analytical processing.

[Nasa MPLnet](https://mplnet.gsfc.nasa.gov/)

##### Afterpulse & Darkcounts

MPLNet Afterpulse & Darkcounts instruction document has been completed and uploaded along with the latex source file in the `mplnet_site_info` directory. Review the [pdf](https://github.com/jdgreene/mplnet-data/blob/main/mplnet_site_info/Afterpulse_Darkcounts.pdf) and update as needed. Image resources for document are also contained within the `mplnet_site_info` folder.

### Requirements
Included `requirements.txt` to assist in setting up Python environment for MPLNet Data Tool
The `mplnet.dat` file is required to be in the same directory as the python application for access to the MPLNet website.

#### Future Work

### Installation Notes

##### Download Python
Python is required and can be downloaded at [Python.org](https://www.python.org/downloads/).
The current version has been tested with Python 3.10.2, 3.11.5, and 3.12.0

##### Environment Setup (optional)
Windows instructions:
Once installed an environment can be setup. Open up a command prompt at the location.
Use the python command `python -m venv env_folder`
Activate the environment from the `Scripts` folder within the `env_folder` : `env_folder/Scripts/activate.bat`
If an environment is not setup the dependencies will be installed system wide.

##### Clone the repo
Clone or download from [github.com](https://github.com/jdgreene/mplnet-data).

##### Install the requirements for the GUI
Using the `requirements.txt` file located in the main folder of the project.
Run the command `pip install -r requirments.txt`

#### Run `mplnet_app.py`
From the command line at the location of the file `mplnet_app.py` run the command `python mplnet_app.py` to start the GUI.
With windows it is often easier to use the file window to navigate to the location of the file then using the mouse right click context menu to open a command window at the location.

