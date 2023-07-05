import netCDF4 as nc
import matplotlib.pyplot as plt 
import numpy as np

# determine the input file path
# file = 'Data/erainterim_moda_L60_N24_200001010000000.nc'
file = 'Data/vsf.data.nc'
# open the file
dataset = nc.Dataset(file)
print(dataset)

