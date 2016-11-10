"""
This modules provides functions that help to:
    * create a single shapefile providing all input data needed to run shyft (model-stack dependend) from
        ** a shapefile holding the target polygons for shyft-cells
        ** a shapefile per basin-unit holding a polygon marking the boundaries of the basin-unit
        ** raster files holding the cell properties
        ** a setup file in yaml format
    * create a shyft-conform inputfile that works in line with the netcdf-repository
"""
# by Felix Matt, University of Oslo, 2016-11-10.

from shapely.geometry import mapping, shape
import fiona
from rasterstats import zonal_stats
import yaml
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal,osr

# intersection algorithm from:
# http://gis.stackexchange.com/questions/74858/fiona-preffered-method-for-defining-a-schema

def get_settings(filename):
    """
    Parameters
    ----------
    filename: str type        
    """
    with open(filename) as stream:
        return yaml.load(stream)

def check_geotiff_crs(setting_file, crs):
    """
    Parameters
    ----------
    setting_file: str type
        yaml settings file. Needed for units etc to map additinal shapefile attributes 
        (e.g. land cover fractions; limited to 10 characters) to longnames in target dataset
    crs: str type
        epsg code of projection; eg.: epsg:32643
    """
    settings = get_settings(setting_file)
    for name,item in settings['properties'].items():
        gtiff = gdal.Open(item['file']) 
        prj = gtiff.GetProjection()
        srs = osr.SpatialReference(wkt=prj)
        crs_tmp = 'epsg:'+srs.GetAttrValue('authority',1)
        if not crs == crs_tmp:
            raise Exception("File {} has different projection than {}: {}".format(item['file'], crs, crs_tmp))

def create_shyft_shapefile(setting_file, outfile):
    """
    Parameters
    ----------
    setting_file: str type
        yaml settings file. Needed for units etc to map additinal shapefile attributes 
        (e.g. land cover fractions; limited to 10 characters) to longnames in target dataset
    outfile: str type
        filename or filename.shp, target filename to write dataset;
        when given without file ending: creates directory and writes files in it with correct file ending
    """
    # load  yaml settings
    settings = get_settings(setting_file)
    # schema of the new shapefile
    schema =  {'geometry': 'Polygon','properties': {'cid':'int', 'area':'float:15.5', 'x':'float:15.5', 'y':'float:15.5'}}
    for name,item in settings['properties'].items():
        schema['properties'][name] = item['data_type']
    crs_dict = {'init': settings['crs']}
    print("schema for writing output: ", schema) 
    # creation of the new shapefile with the intersection
    with fiona.open(outfile, 'w', driver='ESRI Shapefile', schema=schema, crs=crs_dict) as trg:
        for basin_name, basin_dict in settings['basin-units'].items():
            print("Calculating for basin {}...".format(basin_name))
            with fiona.open(basin_dict['file']) as basin_shp:
                crs = basin_shp.crs['init']
                check_geotiff_crs(setting_file, crs)
                for bas in basin_shp:
                    with fiona.open(settings['grid']['file']) as grid_shp:
                        crs_tmp = grid_shp.crs['init']
                        if not crs == crs_tmp:
                            raise Exception("{} and {} have different projections: {} != {}".format(basin_dict['file'], 
                                settings['grid']['file'], crs, crs_tmp))
                        for grid in grid_shp:
                            if shape(bas['geometry']).intersects(shape(grid['geometry'])):
                                shp = shape(bas['geometry']).intersection(shape(grid['geometry']))
                                geometry = mapping(shp)
                                area = shp.area
                                centroid = shp.centroid
                                prop = {'cid':basin_dict['cid'], 'area':area, 'x':centroid.x, 'y':centroid.y}
                                for name,item in settings['properties'].items():
                                    stat = zonal_stats(geometry, item['file'], stats="mean", nodata=settings['nodata'], all_touched=True) # all touched must be set, otherwise returns None if no raster center in polygon.
                                    if not ((stat[0]['mean']>=item['min']) and (stat[0]['mean']<=item['max']) and len(stat)==1):
                                        err_msg = "Resampled value of {} not in range: {} not in range {} - {}"
                                        raise Exception( err_msg.format(name,stat[0]['mean'],item['min'],item['max']) )
                                    prop[name] = stat[0]['mean']
                                trg.write({'geometry':geometry, 'properties': prop})
    print("Output written to {}".format(outfile))

def write_shyft_conform_nc(shapefile, setting_file, outfile):
    """
    Parameters
    ----------
    shapefile: str type
        filename.shp, shapefile filename
        containing cell polygons and attributes (minimum: x,y,z,area, cid)
    setting_file: str type
        yaml settings file. Needed for units etc to map additinal shapefile attributes 
        (e.g. land cover fractions; limited to 10 characters) to longnames in target dataset
    outfile: str type
        filename.nc, target filename to write dataset
    """
    variable_mapping = get_settings(setting_file)['properties']
    with fiona.open(shapefile) as shp:
        nr_cells = len(shp)
        var_dict = {}
        with Dataset(outfile, 'w', format='NETCDF4') as dset:
            cell = dset.createDimension('cell', nr_cells) # only dimension in nr of cells
            
            # every cell needs x, y, area, and catchment_id, and the data needs to be in a certain projection.
            var_dict['x'] = dset.createVariable('x','f8',('cell',))
            var_dict['x'].units = 'm'
            var_dict['x'].axis = 'X'
            var_dict['x'].standard_name = 'projection_x_coordinate'

            var_dict['y'] = dset.createVariable('y','f8',('cell',))
            var_dict['y'].units = 'm'
            var_dict['y'].axis = 'Y'
            var_dict['y'].standard_name = 'projection_y_coordinate'

            var_dict['area'] = dset.createVariable('area','f8',('cell',))
            var_dict['area'].units = 'm^2'
            var_dict['area'].grid_mapping = 'crs'
            var_dict['area'].coordinates = 'y x z'
            
            var_dict['cid'] = dset.createVariable('catchment_id','i4',('cell',))
            var_dict['cid'].units = '-'
            var_dict['cid'].grid_mapping = 'crs'
            var_dict['cid'].coordinates = 'y x z'
            
            # additional variables
            #return variable_mapping 
            for name,item in variable_mapping.items():
                var_dict[name] = dset.createVariable(item['long_name'],'f8',('cell',))
                var_dict[name].unit = item['unit']
                var_dict[name].grid_mapping = 'crs'
                var_dict[name].coordinates = 'y x z'
            
            # fetch values from shapefile fields
            for i,s in enumerate(shp):
                for var_name, var in var_dict.items():
                    var[i] = s['properties'][var_name]
 
            # write crs
            crs_dim = dset.createDimension('crs_dim', None) # only dimension in nr of cells
            crs_var = dset.createVariable('crs','i4',('crs_dim',))
            crs_var.epsg_code = "epsg:{}".format(shp.crs['init'])
            #crs_var.grid_mapping_name = catchment.grid_mapping_name #needed?
            #crs_var.proj4 = catchment.proj4 

if __name__=="__main__":
    setting_file = "example/settings_example.yaml"
    shapefile = "cell-data"
    ncfile = "cell-data.nc"
    create_shyft_shapefile(setting_file, shapefile)
    write_shyft_conform_nc(shapefile, setting_file, ncfile)
