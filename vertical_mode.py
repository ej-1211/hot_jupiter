import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Read in the data
filepath_freq = '/Users/changyichieh/Documents/Hot Jupiter/MODES/Hough_coeff_200001010000000.nc'
data = nc.Dataset(filepath_freq)
print(data)