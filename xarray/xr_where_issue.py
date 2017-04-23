""" xr.where issue: copying a mask from a NetCDF file """
# Copy this file by doing ...
# The NetCDF files are stored in
# ...
import xarray as xr
# Read wind file
f = xr.open_dataset('CCSM4_ens1_19821201_19831130_ws10_0_NAtl_DJFmean.nc')
lat = f.lat
lon = f.lon
ws10 = f.ws10
print(ws10) # To see data
# Read wave file
f = xr.open_dataset('www.Hs.mask.nc')
hs = f.hs
hs = hs.rename({'latitude': 'lat', 'longitude': 'lon'})
hs.coords['lon'] = lon
print(hs) # To see nans
# Try and set nans in the wind file
ws10_masked = ws10.where(hs.isnull())
print(ws10_masked)

