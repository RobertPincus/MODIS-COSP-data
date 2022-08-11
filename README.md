This repository contains Python scripts for manipulating the 
(monthly data)[https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MCD06COSP_M3_MODIS] produced 
by the MODIS cloud group in support of the MODIS simulator for climate models. The data are described in  
Pincus et al, 2022: Updated observations of clouds by MODIS for global model assessment, available at 
(https://crew.ldeo.columbia.edu/people/robert-pincus/)

The data record available at the time of this writing extends from July 2002 to May 2022. 

# Downloading data

The `download-data-with-wget.sh` script attempts to download the entire 16Gb data record to the 
directory specified by the environment variable `MODIS_DATA_CACHE_DIR`. Downloading the data
requires an authentication token for the website in environment variable `EARTHDATA_TOKEN`. 
The simple process for obtaining a token is described on the 
[data distribution website](https://ladsweb.modaps.eosdis.nasa.gov/learn/download-files-using-laads-daac-tokens/). 

At this writing the NASA http server is not always responsive and experiences "Internal Server Error"s a 
small fraction of the time, so it may take several hours for data to be downloaded and the script might need 
to be run more than once (previously downloaded files will be skipped by wget). 

# Post-processing data

As described in the documentation paper, the original form of the data contains variables which 
are not likely to be used by most analysts. The Python script `process-local-data` process data 
for a range of dates (by default, the entire data record) to produce time series of scalar variables 
and properly-normalized joint histograms, which can be written to netCDF files and/or a Zarr store 
(both are enabled by default). The original data are taken from `MODIS_DATA_CACHE_DIR`; the revised 
files are written to the local directory. 