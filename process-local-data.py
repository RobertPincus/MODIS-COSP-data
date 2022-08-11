#! /usr/bin/env python
import os
import datetime
import re
from   pathlib import Path
import xarray as xr
import pandas as pd
import numpy  as np

#
# Post-process the MODIS-COSP dataset described in 
#     Pincus et al, 2022: Updated observations of clouds by MODIS for global model assessment
#    (https://crew.ldeo.columbia.edu/people/robert-pincus/)
#
# The original data come in netcdf4 files with one group per primary variables 
#   This script extracts the space/time mean, per 1-degree grid cell for each month, 
#    and normalizes joint histograms so they are expressed as cloud fractions rather 
#    than a number of observations. 
# The more compact dataset can be written as a single Zarr store and/or as a set 
#    of netCDF files. 
#
write_netcdf = True
write_zarr   = True 

dates = pd.date_range("2002-07-01", "2022-05-01", freq="MS")  # "MS" for "month start"
dates = pd.date_range("2021-01-01", "2021-12-01", freq="MS")  # "MS" for "month start"
cacheDir = Path(os.environ['MODIS_DATA_CACHE_DIR'])

#
# See Table 1 
#
scalar_vars = [
    'Solar_Zenith',
    'Solar_Azimuth',
    'Sensor_Zenith',
    'Sensor_Azimuth',
    'Cloud_Top_Pressure',
    'Cloud_Mask_Fraction',
    'Cloud_Mask_Fraction_Low',
    'Cloud_Mask_Fraction_Mid',
    'Cloud_Mask_Fraction_High',
    'Cloud_Optical_Thickness_Liquid',
    'Cloud_Optical_Thickness_Ice',
    'Cloud_Optical_Thickness_Total',
    'Cloud_Optical_Thickness_PCL_Liquid',
    'Cloud_Optical_Thickness_PCL_Ice',
    'Cloud_Optical_Thickness_PCL_Total',
    'Cloud_Optical_Thickness_Log10_Liquid',
    'Cloud_Optical_Thickness_Log10_Ice',
    'Cloud_Optical_Thickness_Log10_Total',
    'Cloud_Particle_Size_Liquid',
    'Cloud_Particle_Size_Ice',
    'Cloud_Particle_Size_PCL_Liquid',
    'Cloud_Particle_Size_PCL_Ice',
    'Cloud_Water_Path_Liquid',
    'Cloud_Water_Path_Ice',
    'Cloud_Water_Path_PCL_Liquid',
    'Cloud_Water_Path_PCL_Ice',
    'Cloud_Retrieval_Fraction_Liquid',
    'Cloud_Retrieval_Fraction_Ice',
    'Cloud_Retrieval_Fraction_Total',
    'Cloud_Retrieval_Fraction_PCL_Liquid',
    'Cloud_Retrieval_Fraction_PCL_Ice',
    'Cloud_Retrieval_Fraction_PCL_Total',
]


###########################
# Date/filename mapping

def make_filename(date):
    """Make a NetCDF4 filename NASA MODIS-COSP data based on an input date, based on
    the contents of the local directory specified by ``cacheDir``

    :param date: A member of the ``pandas.core.indexes.datetimes.DatetimeIndex``
    """
    day_of_year = date.timetuple().tm_yday
    return ([d for d in cacheDir.glob(f"MCD06COSP_M3_MODIS.A{date.year}{day_of_year:03d}.062.*.nc")][0])
    

###########################
# Scalar variable functions 

def get_from_group(f, g):
    return xr.open_dataset(f, engine="netcdf4", group=g)

def make_time_series(dates, files, vars):
    #
    # Returns a data set with time-mean values of a set of variables
    #
    ds = xr.Dataset(data_vars =
      {v:xr.concat([get_from_group(f, v).Mean.T.assign_attrs(get_from_group(f, v).attrs)
      for f in files], dim="date") for v in vars})
    # Need to delete _FillValue attribute or conflicts will arise when writing 
    for v in vars: 
        for a in ["_FillValue", "scale_factor", "add_offset"]: 
          if a in ds[v].attrs: ds[v].attrs.pop(a)
    #
    # Add coordinates
    #    When accessing a group the lat and lon variables are indexes, not numerical values
    #
    template = xr.open_dataset(files[0], engine="netcdf4")
    ds["latitude"]  = template.latitude
    ds["longitude"] = template.longitude
    ds["date"]      = dates
    # Copy global attributes
    return(ds.assign_attrs(template.attrs))

###########################
# Functions for time series of each joint histogram
def rename_dims(name):
    dimname = re.sub("_[67]",   "", name)    # Remove _6 or _7
    dimname = re.sub("jhisto_", "", dimname) # jhisto_
    #
    # Optical thickness and cloud top pressure bins are the same across joint histograms
    #
    if "optical_thickness" in dimname :
        dimname = re.sub("_total", "", dimname)
        dimname = re.sub("_liquid", "", dimname)
        dimname = re.sub("_ice", "", dimname)
    return(dimname)

def var_name(v, jhisto_name):
    return(v + "_vs_" + jhisto_name)

def ds_name(jhisto_name):
    return("JHisto_vs_" + jhisto_name)

def make_jhisto_series(dates, files, host_var, jhisto):
    #
    # Returns a time series of normalized joint histograms (expressed as cloud fractions)
    #
    group_name = host_var
    if "PCL"    in jhisto: group_name = group_name + "_PCL"
    if "Liquid" in jhisto: group_name = group_name + "_Liquid"
    if "Ice"    in jhisto: group_name = group_name + "_Ice"
    ds = xr.Dataset(data_vars =
      {var_name(host_var, jhisto):
      xr.concat([get_from_group(f, group_name)[ds_name(jhisto)].assign_attrs(get_from_group(f, group_name).attrs) /
                 get_from_group(f, "Cloud_Retrieval_Fraction_Total").Pixel_Counts
      for f in files], dim="date")})
    template = xr.open_dataset(files[0], engine="netcdf4")
    ds["latitude"]  = template.latitude
    ds["longitude"] = template.longitude
    ds["date"]      = dates
    #
    # Add histogram coordinates - central values and bounds
    #   Optical thickness or cloud water path
    #
    dimname    = [d for d in ds.dims if "optical_thickness" in d or "water_path" in d][0]
    boundaries = get_from_group(files[0], group_name)[ds_name(jhisto)].JHisto_Bin_Boundaries
    ds[dimname] = boundaries[:-1] + 0.5 * np.diff(boundaries)
    ds[dimname + "_bnds"] = boundaries # Being written as a coordinate, not sure why
    #
    # Particle size or pressure dimension - the joint variables
    #
    dimname    = [d for d in ds.dims if "particle" in d or "pressure" in d][0]
    boundaries = get_from_group(files[0], group_name)[ds_name(jhisto)].JHisto_Bin_Boundaries_Joint_Parameter
    ds[dimname] = boundaries[:-1] + 0.5 * np.diff(boundaries)
    ds[dimname + "_bnds"] = boundaries # Being written as a coordinate, not sure why
    #
    # Rename dims
    return(ds.rename({d:rename_dims(d) for d in ds.dims if 'jhisto' in d}))

###########################
# Code starts here 
###########################
#
# Scalar quantities (time/space means)
#
paths = [cacheDir.joinpath(make_filename(d)) for d in dates]
print ("Processing ", len(paths), "files ")
all_data = make_time_series(dates, paths, scalar_vars)
print ("Assembled scalar data")
if write_netcdf:
    print ("Writing scalar variables to netcdf")
    all_data.to_netcdf("modis-cosp-scalars.nc", encoding={v:{"zlib": True} for v in scalar_vars})

###########################
#
# Joint histograms - see Table 2
#
#
# Joint histograms with cloud top pressure
#
jhisto = "Cloud_Top_Pressure"
for host_var in ["Cloud_Optical_Thickness_" + pcl + phase
                 for pcl in ["", "PCL_"] for phase in ["Total", "Liquid", "Ice"]]:
    hname = f"{host_var}_vs_{jhisto}"
    print (f"Making joint histogram {hname}")
    temp = make_jhisto_series(dates, paths, host_var, jhisto)
    all_data[hname] = temp[hname]
    for v in temp.variables:
        if "_bnds" in v: all_data[v] = temp[v] 
    if write_netcdf: 
        print ("  Writing to netcdf")
        temp.to_netcdf(f"modis-cosp-{hname}.nc", encoding={hname:{"zlib": True}})


#
# Joint histograms with particle size
#
for host_var in ["Cloud_Optical_Thickness", "Cloud_Water_Path"]:
    for pcl in ["", "PCL_"]:
        for phase in ["Liquid", "Ice"]:
            jhisto = f"Cloud_Particle_Size_{pcl}{phase}"
            hname = f"{host_var}_vs_{jhisto}"
            print (f"Making joint histogram {hname}")
            temp = make_jhisto_series(dates, paths, host_var, jhisto)
            all_data[hname] = temp[hname]
            for v in temp.variables:
                if "_bnds" in v: all_data[v] = temp[v] 
            if write_netcdf: 
              print ("  Writing to netcdf")
              temp.to_netcdf(f"modis-cosp-{hname}.nc", encoding={hname:{"zlib": True}})

if write_zarr: all_data.to_zarr("modis-cosp.zarr")
