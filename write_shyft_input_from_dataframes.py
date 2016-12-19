"""
This modules provides the functionality of writing SHyFT conform netCDF data files from provided pandas DataFrames,
holding discharge and forcing data. Information about the specific stations must be provided in a dict-type staion-info file.
"""
# by Felix Matt, University of Oslo, 2016-12-19.

import numpy as np
import calendar as cal
from netCDF4 import Dataset
from os import path

def check_time_naive(datetime_lst):
    """
    checks is all of the datetime objects in the input list are naive, meaning that no timezone information is added.
    If this is the case, we assume the time to be utc, and the python calendar module can be used to calculate timestamps.
    If the time is not naive, we stop the calculation. Time zone handling is not icluded.
    """
    for date in datetime_lst:
        if date.tzinfo is not None:
            raise ValueError("There is information about the timezone. We expect naive UTC datetime objects to calculate the UTC timestamp.")

def write_q_df(dset, station_info, idx_map, name_lst):
    """
    Writing discharge specific variables to existing nc database.
    Parameters:
    -----------
    dset, netCDF4 Dataset type, write access.
    station_info, dict type.
    idx_map, list type, holding idx of items to write to database.
    name_lst, list type, station names fitting to idx_map.
    """
    vlen_t = dset.createVLType(np.int32, "vlen") # variable length type for sequence of arrays with variable lengths
    catchment_id = dset.createVariable('catchment_id', vlen_t, ('series',))
    for i,(idx,name) in enumerate(zip(idx_map,name_lst)):
        catchment_id[i] = station_info['catchment_id'][idx]

def write_met_df(dset, station_info, idx_map, name_lst):
    """
    Writing meteo specific variables to existing nc database.
    Parameters:
    -----------
    dset, netCDF4 Dataset type, write access.
    station_info, dict type.
    idx_map, list type, holding idx of items to write to database.
    name_lst, list type, station names fitting to idx_map.
    """
    z = dset.createVariable('z', np.float64, ('series',))
    z.units = 'm'
    z.axis = 'Z'
    z.standard_name = 'height'
    z.long_name = 'height above mean sea level'
    for i,(idx,name) in enumerate(zip(idx_map,name_lst)):
        z[i] = station_info['z'][idx]

def write_df_to_nc(df, station_info, target_path):
    """
    Write SHyFT conform database from pandas data frames and station info dict.
    Parameters:
    -----------
    df, dict type, keyed with variable name, holding pandas data frame types keyed with station names.
    station_info, dict type, ditionary, keyed with variable name,
    holding dictionry keys with specifics about variable's stations.
    target_path, outfile path (filenames prescribed via variable name).
    """
    
    # some assertions first
    kind = station_info['kind']
    filename = path.join(target_path,kind+'.nc')
    print("Writing {} to SHyFT conform nc file...".format(kind))
    
    # get stations listed in df AND station_info:
    info_name = station_info['name']
    name_lst = []
    idx_map = []
    for key in df.keys():
        if key in info_name:
            name_lst.append(key)
            idx_map.append(info_name.index(key))
    
    # wirte database
    with Dataset(filename, 'w', format='NETCDF4') as dset:
        # dimensions
        series_dim = dset.createDimension('series', size=len(name_lst))
        time_dim = dset.createDimension('time', size=df.index.size)
        
        # variables
        time = dset.createVariable('time', np.float64, ('time',))
        time.units = 'seconds since 1970-01-01 00:00:00 +00:00'

        series_name = dset.createVariable('series_name', str, ('series',))
        series_name.cf_role = 'timeseries_id'

        crs = dset.createVariable('crs', np.int32, ())
        crs.grid_mapping_name = 'transverse_mercator'
        crs.epsg_code = "epsg:{}".format(station_info['epsg'])
        zone = str(station_info['epsg'])[-2:]
        crs.proj4 = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(zone) # Assumes utm on WGS84

        x = dset.createVariable('x', np.float64, ('series',))
        x.units = 'm'
        x.axis = 'X'
        x.standard_name = 'projection_x_coordinate'
        
        y = dset.createVariable('y', np.float64, ('series',))
        y.units = 'm'
        y.axis = 'Y'
        y.standard_name = 'projection_y_coordinate'
        
        ts = dset.createVariable(kind, np.float64, ('time', 'series'))
        ts.units = station_info['units']
        ts.grid_mapping = 'crs'
        ts.coordinates = 'y x z'

        # fill variables
        check_time_naive(df.index)
        utc_time = [cal.timegm(date.timetuple()) for date in df.index] # calculates utc time stamp from utc datetime object
        time[:] = utc_time
        series_name[:] = np.array(name_lst)
        
        for i,(idx,name) in enumerate(zip(idx_map,name_lst)):
            x[i] = station_info['x'][idx]
            y[i] = station_info['y'][idx]
            ts[:,i] = df[name].values
        
        # "kind" specific variables
        if kind == 'discharge':
            write_q_df(dset, station_info, idx_map, name_lst)
        else:
            write_met_df(dset, station_info, idx_map, name_lst)
    print("... done.")
