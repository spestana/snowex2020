"""
Functions for working with Airborne TIR imagery from SnowEx 2020 Grand Mesa.
Steven Pestana (spestana@uw.edu)
"""

import scipy.io as sio
import mat73 # https://github.com/skjerns/mat7.3
import xarray as xr
import rioxarray
import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio as rio
import rasterio.plot as rioplt
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling

def open_parse_geotiffs(files):
    '''Given a list of filepaths to GeoTiff images, opens, parses timestamps from the filenames, appends to list'''
    timestamps_UTC = []
    dataarrays = []
    
    for file in files:
        # get timestamp from filename
        datetimelist = file.split("\\")[-1].split("_")[-1].split(".")[0].replace("T","-").split("-")
        #print('parsing: ', datetimelist)
        y = int(datetimelist[0])
        mo = int(datetimelist[1])
        d = int(datetimelist[2])
        h = int(datetimelist[3][:2])
        m = int(datetimelist[3][2:4])
        s = int(datetimelist[3][4:])
        timestamps_UTC.append(pd.Timestamp(y, mo, d, h, m, s))
        #print(airborne_ir_timestamp_UTC)
        
        # open file
        dataarrays.append(rioxarray.open_rasterio(file))
        
    return timestamps_UTC, dataarrays

def ir_mat2dataset(mat_filepath):
    ''' Given a filepath to an airborne IR mosaic .mat file, converts to an xarray dataset.'''

    # read in .mat file
    aircraft_data_mat = sio.loadmat(mat_filepath)  
    
    # create an xarray dataset from the data
    # create dataset
    ds = xr.Dataset({
        'STCtemp': xr.DataArray(
                    data   = aircraft_data_mat['STCtemp'], # scaled IR temperature data
                    dims   = ['y', 'x', 'time'],
                    coords = {'x': aircraft_data_mat['Eutm'].squeeze(), # easting coordinates (MGRS-UTM, m
                              'y': aircraft_data_mat['Nutm'].squeeze(), # northing grid coordinates (MGRS-UTM, 
                              'time': aircraft_data_mat['time'].squeeze()}, # time in epoch
                    attrs  = {
                        '_FillValue'  : np.nan,
                        'units'       : 'C',
                        'description' : 'Airborne thermal infrared from University of Washington Applied Physics Lab Compact Airborne System for Imaging the Environment instrument.'
                        }
                    ),
        'zDEM': xr.DataArray(
                    data   = aircraft_data_mat['zDEM'], # elevation (meters)
                    dims   = ['y', 'x'],
                    coords = {'x': aircraft_data_mat['Eutm'].squeeze(),
                              'y': aircraft_data_mat['Nutm'].squeeze()},
                    attrs  = {
                        '_FillValue'  : np.nan,
                        'units'       : 'm',
                        'description' : 'DEM from SRTM 1 arc-second product merged, and clipped around Grand Mesa https://doi.org/10.5066/F7PR7TFT'
                        }
                )
    },
            attrs = {'example_attr': 'this is a global attribute'}
        )
    
   
    # convert matlab time to python datetime (UTC)
    M = ds.time.values
    base_date = datetime.datetime(1970,1,1,0,0,0)
    return_matrix = []
    for time in M:
        return_matrix.append(base_date + datetime.timedelta(seconds=time))
    ds['time'] = return_matrix
    
    # set coordinate system
    ds.rio.set_spatial_dims('y', 'x', inplace=True)
    ds.rio.write_crs("epsg:32613", inplace=True) # UTM Zone 13N
    
    return ds


def eo_mat2dataset(mat_filepath):
    ''' Given a filepath to an airborne EO (visible imagery) mosaic .mat file, converts to an xarray dataset.'''
    
    aircraft_eo_data_mat = mat73.loadmat(mat_filepath)
    
    # convert time in epoch to isoformat
    #time = [datetime.datetime.utcfromtimestamp(this_time).isoformat() for this_time in aircraft_eo_data_mat['time'].squeeze()]

    
    # scale/stretch RGB colors to values between 0 and 1
    original_max = np.nanmax(aircraft_eo_data_mat['SRGB'])
    original_min = np.nanmin(aircraft_eo_data_mat['SRGB'])
    scaled_RGB_data = (aircraft_eo_data_mat['SRGB'] - original_min) / (original_max - original_min)
    #scaled_RGB_data = scaled_RGB_data.astype('int')
    scaled_RGB_data[np.isnan(scaled_RGB_data)] = 0
    
    # create dataset
    ds = xr.Dataset({
        'SRGB': xr.DataArray(
                    data   = scaled_RGB_data, # Red, Green, Blue image bands
                    dims   = ['y', 'x', 'band', 'time'],
                    coords = {'x': aircraft_eo_data_mat['Xs'][0,:], # easting coordinates (MGRS-UTM, m)
                              'y': aircraft_eo_data_mat['Ys'][:,0], # northing grid coordinates (MGRS-UTM, m)
                              'band': ['r', 'g', 'b'], # labels for red, green, blue image bands
                              'time': aircraft_eo_data_mat['time'].squeeze()}, # time in isoformat
    
                    attrs  = {
                        '_FillValue'  : np.nan,
                        'units'       : 'digital numbers',
                        'description' : 'visible imagery',
                        'timezone'    : 'time in UTC'
                        }
                    ),
        'zDEM': xr.DataArray(
                    data   = aircraft_eo_data_mat['zDEM'], # elevation (meters)
                    dims   = ['y', 'x'],
                    coords = {'x': aircraft_eo_data_mat['Xs'][0,:], # easting coordinates (MGRS-UTM, m)
                              'y': aircraft_eo_data_mat['Ys'][:,0]}, # northing grid coordinates (MGRS-UTM, m)
                    attrs  = {
                        '_FillValue'  : np.nan,
                        'units'       : 'meters',
                        'description' : 'digital elevation map used for georectification',
                        'source'      : 'https://doi.org/10.5066/F7PR7TFT',
                        'horizontal datum' : 'WGS84',
                        'vertical datum' : 'EGM96',
                        'original resolution' : '30 m / 1 arc-second'
                        }
                )
    },
            attrs = {'description': 'Airborne visible imagery from University of Washington Applied Physics Lab Compact Airborne System for Imaging the Environment instrument. Imagery of Grand Mesa, Colorado from SnowEx 2020 field campaign.'}
        )
    
    
    # convert matlab time to python datetime (UTC)
    M = ds.time.values
    base_date = datetime.datetime(1970,1,1,0,0,0)
    return_matrix = []
    for time in M:
        return_matrix.append(base_date + datetime.timedelta(seconds=time))
    ds['time'] = return_matrix
    
    # set coordinate system
    ds.rio.set_spatial_dims('y', 'x', inplace=True)
    ds.rio.write_crs("epsg:32613", inplace=True) # UTM Zone 13N / NAD83
    
    return ds

def ir_zonal_stats(ir_filepath, shapefile_filepath, shapefile_number, return_masked_array=False):
    '''Calculate zonal statistics for an ASTER TIR geotiff image within a single polygon from a shapefile.'''
    with rio.open(ir_filepath) as src:
        
        # Open the shapefile
        zone_shape = gpd.read_file(shapefile_filepath)

        # Make sure our shapefile is the same CRS as the IR image
        zone_shape = zone_shape.to_crs(src.crs)
        
        # Mask the IR image to the area of the shapefile
        try:
            masked_ir, mask_transform = mask(dataset=src, 
                                            shapes=zone_shape.geometry[shapefile_number],
                                            crop=True,
                                            all_touched=True,
                                            filled=True)
        # Note that we still have a "bands" axis (of size 1) even though there's only one band, we can remove it below
        except ValueError as e: 
            # ValueError when shape doesn't overlap raster
            print(e)
            return
        
        # change data type to float64 so we can fill in any values =0 with NaN values
        masked_ir = masked_ir.astype('float64')
        masked_ir[masked_ir==0] = np.nan
                
        
        # Remove the extra dimension (bands, we only have one band here)
        masked_ir_tb = masked_ir.squeeze()
        # Get all pixel values in our masked area
        values = masked_ir_tb[0,:,:].flatten() # flatten to 1-D
        values = values[~np.isnan(values)] # remove NaN pixel values
        
        # Calculate zonal statistics for this area (mean, max, min, std:)
        try:
            masked_ir_tb_mean = values.mean()
            masked_ir_tb_max = values.max()
            masked_ir_tb_min = values.min()
            masked_ir_tb_std = values.std()
        except ValueError as e:
            # ValueError when the shapefile is empty I think
            print(e)
            return
        
        if len(values) == 0:
            print('No valid pixels within shapefile bounds')
            return None
        
        if return_masked_array == True:
            return masked_ir_tb_mean, masked_ir_tb_max, masked_ir_tb_min, masked_ir_tb_std, masked_ir_tb
        else:
            return masked_ir_tb_mean, masked_ir_tb_max, masked_ir_tb_min, masked_ir_tb_std