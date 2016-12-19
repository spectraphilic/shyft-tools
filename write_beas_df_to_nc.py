"""
This modules provides the functionality to write data, provided as pandas DataFrames, to SHyFT conform netCDF4 files.
"""
# by Felix Matt, University of Oslo, 2016-12-19.

import numpy as np
import csv
import pandas as pd
import glob

# local
from fetch_variables_from_xls import get_workbook
from fetch_variables_from_xls import fetch_dataframe
from fetch_variables_from_xls import fetch_discharge_ts
from fetch_variables_from_xls import get_workbook 
from fetch_variables_from_xls import get_station_info_from_csv 
from write_shyft_input_from_dataframes import write_df_to_nc
from fetch_radiation_from_WFDEI import get_WFDEI_radiation_for_beas

def get_beas_data():
    """
    Reads from collection of data locations.
    Calls functions to fetch Beas dataframes as variable dicts and stations info dicts keys with variable.
    Some variable specific conversions are done, eg. precip in mm/h and mean temp from max and min temp.
    """
    df_dict = {}
    station_info = {}

    # precip
    filename='/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/rainfall-2011/allrains_1979-2011.xlsx'
    sheets = ['Banjar', 'Bhuntar', 'Larji', 'Pandoh',  'Sainj', 'Manali', 'Janjehli']
    kind = "precipitation"
    df_dict[kind] = fetch_dataframe(filename, sheets=sheets, kind=kind) / 24.
    print("Converted {} from daily sum to mm/h.".format(kind))
    station_info_file = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/stations/precipitation_stations_utm-43n.csv'
    station_info[kind] = get_station_info_from_csv(station_info_file)

    # temp max
    filename = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/temp_humid-2005/Max Temp71-05.xlsx'
    sheets = ['pandoh', 'Bhuntar', 'Larji', 'Manali']
    kind = "temperature"
    df_temp_max = fetch_dataframe(filename, sheets=sheets, kind=kind)

    # temp min
    filename = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/temp_humid-2005/Min Temp71-05.xlsx'
    sheets = ['pandoh', 'Bhuntar', 'Larji', 'Manali']
    kind = "temperature"
    df_temp_min = fetch_dataframe(filename, sheets=sheets, kind=kind)

    # temperature
    df_dict[kind] = (df_temp_min+df_temp_max)/2.
    station_info_file = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/stations/temperature_stations_utm-43n.csv'
    station_info[kind] = get_station_info_from_csv(station_info_file)

    # relative-humidity
    filename = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/temp_humid-2005/humidity71-05.xlsx'
    sheets = ['pandoh', 'Bhuntar', 'Larji', 'Manali']
    kind = "relative_humidity"
    df_dict[kind] = fetch_dataframe(filename, sheets=sheets, kind=kind) / 100.
    print("Converted {} from precent to fraction.".format(kind))
    station_info_file = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/stations/relative_humidity_stations_utm-43n.csv'
    station_info[kind] = get_station_info_from_csv(station_info_file)

    # discharge
    filenames = glob.glob("/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/dischage-2011/*xls")
    names = [fn.split("/")[-1].split(' ')[-1].split('-')[-2] for fn in filenames]
    kind = "discharge"
    df_dict[kind] = fetch_discharge_ts(filenames, names, kind=kind)
    station_info_file = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/stations/discharge_stations_utm-43n.csv'
    station_info[kind] = get_station_info_from_csv(station_info_file)

    # radiation
    kind = 'global_radiation'
    df_dict[kind], station_info[kind] = get_WFDEI_radiation_for_beas()

    return station_info, df_dict

def write_beas_data_to_nc(outfile):
    """
    Calls functions to write Beas data frames dict to nc files.
    Parameters:
    -----------
    outfile -> location of output file.
    """
    station_info, df_dict = get_beas_data()
    for key in station_info.keys():
        write_df_to_nc(df_dict[key], station_info[key], outfile)

if __name__=='__main__':
    """Fetches Beas data from source and writes it to SHyFT conform netCDF4 files."""
    write_beas_data_to_nc("../shyft_input_nc/")

