"""
This modules provides functionality to read data from the WFDEI database (Watch Forcing Data Era Interim)
as a pandas DataFrame. The approximate region of which data is extraced is specified by a Bounding Box in lat and lon.
"""
# by Felix Matt, University of Oslo, 2016-12-19.

import numpy as np
import glob
from netCDF4 import Dataset
import datetime
import pandas as pd
import pyproj

def get_bounding_idxs(lats, lons, bounding_box):
    """
    Returns boundary indexes for lat, lon inside bounding box.
    parameters:
    -----------
    lat, np.array type, 1-d vector.
    lon, np.array type, 1-d vector.
    bounding_box, list type, example = [lat_min, lat_max, lon_min, lon_max].
    """
    idxs = []
    idxs.append(np.argmin(abs(lats-bounding_box[0]))-1) # be bit bigger
    idxs.append(np.argmin(abs(lats-bounding_box[1]))+2) # be bit bigger + array idx selection in python
    idxs.append(np.argmin(abs(lons-bounding_box[2]))-1)
    idxs.append(np.argmin(abs(lons-bounding_box[3]))+2)
    return idxs

def get_timestamp(time):
    """
    parameters:
    -----------
    times, WFDEI 'time' variable type
    return datetime vector in UTC
    """
    start_date = datetime.datetime.strptime(time.units[14:], '%Y-%m-%d %H:%M:%S')
    dates = [start_date+datetime.timedelta(seconds=float(sec)) for sec in time[:]]
    return np.array(dates)

def get_WFDEI_elevations(filename, bounding_idx):
    """provides WFDEI elevation data from netCDF databasei (np.array type)."""
    with Dataset(filename) as dset:
        elevation = dset.variables['elevation'][bounding_idx[0]:bounding_idx[1],
                                                bounding_idx[2]:bounding_idx[3]]
    return elevation

def get_WFDEI_variable(filenames, varname, bounding_idx):
    """
    Fetch variable from WFDEI database, region on globe specified with Bounding Box.
    parameters:
    -----------
    filenames, list type, sorted list with netcdf filenames
    varname, str type, WFDEI variable name
    bounding_idx, list type, bounding_box indicees in global variable array,
    e.g.: [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max].
    """
    print("Reading variable {} ...".format(varname))
    for i,fn in enumerate(filenames):
        with Dataset(fn) as dset:
            var = dset.variables[varname][:,bounding_idx[0]:bounding_idx[1],
                                            bounding_idx[2]:bounding_idx[3]]
            t = dset.variables['time']
            if i==0:
                arr = var + 0. # make copies
                time = get_timestamp(t) + datetime.timedelta(0)
            else:
                arr = np.concatenate((arr,var), axis=0)
                time = np.concatenate((time,get_timestamp(t)))
        print("{}% .. ".format(int(float(i+1)/len(filenames)*100)), end='')
    print(" DONE!")
    return time, arr

def make_data_frame_and_station_info(dates, var, elevation, lat, lon, ref_system, epsg, kind, units):
    """
    Converts data array in pandas DataFrame and station info
    which can be used to write to SHyFT conform input data.
    Input data is expected to be in lat lon projection, target projection can be specified via epsg code.
    parameters:
    -----------
    dates, array type (vector).
    var, array type (3d).
    elevation, array type (2d).
    lat, array type (vector), latitude of grid points.
    lon, array type (vector), longitude of grid points.
    ref_system, str type, e.g. "utm-43n".
    epsg, int type, epsg code of target projection.
    kind, str type, variable name (e.g. "radiation").
    units, str type, unit specification (e.g. "m s-1").
    returns pandas DataFrame type of station data (keyed with station name).
    returns dict type of station info of each station in DataFrame.
    """
    proj = pyproj.Proj(init='epsg:'+str(epsg))
    df = {}
    station_info = {}
    station_info['kind'] = kind
    station_info['ref_system'] = ref_system
    station_info['epsg'] = epsg
    station_info['units'] = units
    station_info['name'],station_info['x'],station_info['y'],station_info['z'] = [],[],[],[]
    for i,la in enumerate(lat):
        for j,lo in enumerate(lon):
            name = 'station_'+str(i)+'_'+str(j)
            df[name] = pd.Series(var[:,i,j], index=dates)
            x, y = proj(lo,la)
            station_info['name'].append(name)
            station_info['x'].append(x)
            station_info['y'].append(y)
            station_info['z'].append(elevation[i,j])

    return pd.DataFrame(df), station_info

def get_WFDEI_radiation_for_beas():
    """ Fetches pandas DataFrame and station info with Beas catchment specific bounding box and projection (utm-43n)."""
    files = sorted(glob.glob("/home/felixm/hycamp/SHYFT_DATA/Global/WFDEI_3hrs/data/*.nc"))
    bounding_box = [31.3, 32.5, 76.8, 78.0]
    #bounding_box = [31.3, 32.5, 74.8, 80.0] # test
    with Dataset(files[0]) as dset:
        lat,lon = dset.variables['lat'][:], dset.variables['lon'][:]
    bounding_idx = get_bounding_idxs(lat, lon, bounding_box)
    elevations = get_WFDEI_elevations("/home/felixm/hycamp/SHYFT_DATA/Global/WFDEI_3hrs/WFDEI-elevation.nc", bounding_idx)
    dates, var = get_WFDEI_variable(files, "SWdown", bounding_idx)
    lat_sel = lat[bounding_idx[0]:bounding_idx[1]]
    lon_sel = lon[bounding_idx[2]:bounding_idx[3]]
    epsg = 32643
    ref_system = 'utm-43n'
    df, station_info = make_data_frame_and_station_info(dates, var, elevations, 
            lat_sel, lon_sel, ref_system, epsg, "global_radiation", "W m-2")
    return df, station_info

if __name__ == '__main__':
    df, station_info = get_WFDEI_radiation_for_beas()
