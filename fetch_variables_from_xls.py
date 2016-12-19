"""
This modules provides the functionality of reading data from a xls database. The readers are specified to read from the Beas catchment database. Additinal station information is read from csv tables.
"""
import numpy as np
from openpyxl import load_workbook
import calendar as cal
import pandas as pd
from IPython import embed

def get_station_info_from_csv(filename):
    """Convert station information from csv file to dict (return)."""
    station_info = {}
    station_info_raw = np.genfromtxt(filename, names=True, delimiter=',', dtype=None)
    # some assertions first
    assert(len(set(station_info_raw['kind']))==1), "Different kinds of variables in station_info."
    assert(len(set(station_info_raw['ref_system']))==1), "Different kinds of reference systems in station_info."
    assert(len(set(station_info_raw['epsg']))==1), "Different kinds of epsg codes in station_info."
    assert(len(set(station_info_raw['units']))==1), "Different kinds of units in station_info."

    station_info['kind'] = station_info_raw['kind'][0].decode('UTF-8')
    station_info['ref_system'] = station_info_raw['ref_system'][0].decode('UTF-8')
    station_info['epsg'] = int(station_info_raw['epsg'][0])
    station_info['units'] = station_info_raw['units'][0].decode('UTF-8')
    station_info['name'] = [name.decode('UTF-8') for name in station_info_raw['name']]
    station_info['x'] = [int(xval) for xval in station_info_raw['x']]
    station_info['y'] = [int(yval) for yval in station_info_raw['y']]
    if 'z' in station_info_raw.dtype.names:
        station_info['z'] = [int(yval) for yval in station_info_raw['z']]
    if 'catchment_id' in station_info_raw.dtype.names:
        catchment_id = []
        for catch_id in station_info_raw['catchment_id']:
            catchment_id.append(np.array([np.int32(cid) for cid in catch_id.decode('UTF-8').split('-')]))
        station_info['catchment_id'] = catchment_id
    return station_info

def get_workbook(filename):
    """Open xls workooks."""
    return load_workbook(filename=filename, read_only=True)

def fetch_precip_ts(df):
    """
    Postprocess dataframe of precipitation data read from
    xls file specific to the Beas data excel layout.
    Parameters:
    -----------
    df, pandas dataframe type.
    """
    df = df.replace(0,-9999)
    df = df.replace(np.nan,0)
    df = df.replace(-9999,np.nan)
    all_vals = []
    years = np.empty(df.shape[1])
    pops_noleap = [380, 348, 317, 285, 254, 222, 190, 159, 127, 96, 64, 63, 34, 2] # lines in excel sheet to pop
    pops_leap = [380, 348, 317, 285, 254, 222, 190, 159, 127, 96, 64, 34, 2]
    for i,(y,col) in enumerate(df.items()):
        year = 1900+y if y<2000 else y
        #col_lst = list(col[:379])
        col_lst = [float(val) if type(val) is not str else np.nan for val in col[:379]]
        if cal.isleap(year):
            for idx in pops_leap:
                col_lst.pop(idx-2)
        else:
            for idx in pops_noleap:
                col_lst.pop(idx-2)
        all_vals += col_lst
        years[i] = year
        if year == 2011: # don't go further than 2011 column
            break
    dates = pd.date_range(str(int(round(years[0])))+'-01-01',str(int(round(years[-1])))+'-12-31')
    return pd.Series(all_vals, index=dates)
    
def fetch_temperature_ts(df):
    """
    Postprocess dataframe of temperature data read from 
    xls file specific to the Beas data excel layout.
    Parameters:
    -----------
    df, pandas dataframe type.
    """
    all_vals = []
    years = np.empty(df.shape[1])
    pops_noleap = [62,2] # lines in excel sheet to pop
    pops_leap = [2]
    for i,(y,col) in enumerate(df.items()):
        year = 1900+y if y<2000 else y
        col_lst = [float(val) if type(val) is not str else np.nan for val in col]
        if cal.isleap(year):
            for idx in pops_leap:
                col_lst.pop(idx-2)
        else:
            for idx in pops_noleap:
                col_lst.pop(idx-2)
        all_vals += col_lst
        years[i] = year
    dates = pd.date_range(str(int(round(years[0])))+'-01-01',str(int(round(years[-1])))+'-12-31')
    return pd.Series(all_vals, index=dates)

def fetch_discharge_ts(filenames, names, kind='discharge'):
    """
    Read and postprocess dataframe of discharge data read from 
    xls file specific to the Beas data excel layout.
    Parameters:
    -----------
    filenames, list type
        list of filenames.
    names, list type
        list of catchments names, same length as filenames.
    kind, str type
        'discharge' or 'gauge'. Deflault = 'discharge'.
    """
    series_dict = {}
    for filename,name in zip(filenames,names):
        df = pd.read_excel(filename, header=1)
        all_vals = []
        assert(df.shape[1]%2==0), "Header does not have two rows er year." 
        years = np.empty(int(round(df.shape[1]/2)))
        pops_noleap = [64,4,3] # lines in excel sheet to pop
        pops_leap = [4,2]
        for i,keys in enumerate(zip(df.keys()[0::2], df.keys()[1::2])):
            if kind == 'discharge':
                col = df[keys[1]]
            elif kind == 'gauge':
                col = df[keys[0]]
            y = keys[0]
            year = 1900+y if y<2000 else y
            col_lst = [float(val) if type(val) is not str else np.nan for val in col]
            if cal.isleap(year):
                for idx in pops_leap:
                    col_lst.pop(idx-2)
            else:
                for idx in pops_noleap:
                    col_lst.pop(idx-2)
            all_vals += col_lst
            years[i] = year
        dates = pd.date_range(str(int(round(years[0])))+'-01-01',str(int(round(years[-1])))+'-12-31')
        series_dict[name] = pd.Series(all_vals, index=dates)
    return pd.DataFrame(series_dict)

def fetch_dataframe(filename, sheets=None, kind=None):
    """
    Fetch dataframe from xls file of Beas meteorologic variable.
    Parameters:
    -----------
    filename, str type, file to read from.
    sheets, list type, holding strings with sheet names.
    kind, str type, meteorologic variable to read.
    """
    print("fetching {} from file {} ...:".format(kind, filename))
    ts_dict = {}
    for sheet in sheets:
        print(sheet)
        df = pd.read_excel(filename, sheet)
        if kind == "precipitation":
            ts_dict[sheet] = fetch_precip_ts(df)
        elif kind == "temperature" or "relative-humidity":
            ts_dict[sheet] = fetch_temperature_ts(df)
        else:
            print("kind not known...")
    return pd.DataFrame(ts_dict)


if __name__=="__main__":
    filename = '/home/felixm/00_PhD/enki/enki_data/beas/MODDRFS/OriginalBeasData/rainfall-2011/allrains_1979-2011.xlsx'
    sheets = ['Banjar', 'Bhuntar', 'Larji', 'Pandoh',  'Sainj', 'Manali']
    df = fetch_dataframe(filename, sheets=sheets, kind='precipitation')

