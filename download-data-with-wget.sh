#! /bin/bash
#
# Download the entire monthly MODIS COSP data plus one example daily file 
#   This shell script is based on the wget script shown at https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/62/MCD06COSP_M3_MODIS/
#   under "See wget script"
#   The data distribution servers require user registration and a token from the 
#   web site following https://ladsweb.modaps.eosdis.nasa.gov/learn/download-files-using-laads-daac-tokens/
#   This script retrieves an authentication token from the environment variable EARTHDATA_TOKEN
#   
# The data are downloaded into a single flat directory specified by environment variable MODIS_DATA_CACHE_DIR
#
# At this writing the NASA http server is not always responsive and experiences "Internal Server Error"s a 
#   small fraction of the time, so it may take several hours for data to be downloaded and the script might need 
#   to be run more than once (previously downloaded files will be skipped by wget)
#
#  The downloaded data may then be post-processed with 'process-local-data.py' to remove unneeded fields and 
#     normalize joint histograms (convert from counts to cloud fractions)
#  See Pincus et al, 2023: Updated observations of clouds by MODIS for global model assessment
#    (https://crew.ldeo.columbia.edu/people/robert-pincus/)
#
# Sample daily file
wget -e robots=off -m -nv -np -R .html,.tmp,.csv,.json -nd --cut-dirs=1 \
   "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/62/MCD06COSP_D3_MODIS/2021/196/" \
   --header "Authorization: Bearer ${EARTHDATA_TOKEN}" -P ${MODIS_DATA_CACHE_DIR}

# All monthly files, year-by-year 
for yr in {2002..2023..1}
do
  echo "Downloading $yr"
  wget -e robots=off -m -nv -np -R .html,.tmp,.csv,.json -nd --cut-dirs=2 \
  "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/62/MCD06COSP_M3_MODIS/${yr}/" \
  --header "Authorization: Bearer ${EARTHDATA_TOKEN}" -P ${MODIS_DATA_CACHE_DIR}
done