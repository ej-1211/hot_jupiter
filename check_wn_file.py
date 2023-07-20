import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from analyzing_hj_mode import process_freq,find_eigenfrequncy

# =======================
# 0. Set which zonal wave number to plot
target_hj = 'WASP 19 b'
zonal_wave_number = 2 #k : 0 means the "first" zonal mode but the index is 0 ; 41 in total (0~40)
vertical_wave_numebr = 1 #n : 1 means the first vertical mode ; 0 means nothing ; 20 in total 
meridional_wave_number = 2 #m : 1 means the first meridional mode ; 0 means nothing ; 30 in total
meridional_mode_type = 'ROT' # 'EIG' or 'WIG' or 'ROT'

# for eigenfrequency:
print_whole_freq_for_all_vertical_modes = True
# =======================
# 1. Read in the data
filepath_wn = f"hot_jupiter/MODES/Data/Output for {target_hj}/hough.nc/hough.wn{zonal_wave_number:05d}.nc"
filepath_freq = 'hot_jupiter/MODES/Data/freq_0718.dat'
wn_data = nc.Dataset(filepath_wn)
freq_data = process_freq(filepath_freq)
# print(wn_data)

# 2. Find the eigenfrequency
if print_whole_freq_for_all_vertical_modes:
    # print the whole eigenfrequency for all vertical modes
    eigenfreq_list = []
    for i in range(20):
        eigenfreq = find_eigenfrequncy(zonal_wave_number, i+1, meridional_mode_type, meridional_wave_number)
        eigenfreq_list.append(eigenfreq)
    print('--------------------------------------------------')
    print(f'The Eigenfrequency of zonal wave number {zonal_wave_number}, vertical wave number 01~20, meridional type {meridional_mode_type} wave number {meridional_wave_number} is')
    for i in [0,5,10,15]:
        print(f'{i+1:02d} : {eigenfreq_list[i]:.10f}\t{i+2:02d} : {eigenfreq_list[i+1]:.10f}\t{i+3:02d} : {eigenfreq_list[i+2]:.10f}\t{i+4:02d} : {eigenfreq_list[i+3]:.10f}\t{i+5:02d} : {eigenfreq_list[i+4]:.10f}')
    print('--------------------------------------------------')
elif print_whole_freq_for_all_vertical_modes == False:
    eigenfreq = find_eigenfrequncy(zonal_wave_number, vertical_wave_numebr, meridional_mode_type, meridional_wave_number)
    print('--------------------------------------------------')
    print(f'The Eigenfrequency of zonal wave number {zonal_wave_number}, vertical wave number {vertical_wave_numebr}, meridional type {meridional_mode_type} wave number {meridional_wave_number} : {eigenfreq}')
    print('--------------------------------------------------')
# 2. plot the uvz profile
plt.figure(figsize=(6, 10))
lat = wn_data['lat'][:]
if meridional_mode_type == 'ROT': meridional_mode_type = 'BAL' # change the name due to the data stucture
u = wn_data[meridional_mode_type][vertical_wave_numebr-1][0].T[meridional_wave_number-1]
v = wn_data[meridional_mode_type][vertical_wave_numebr-1][1].T[meridional_wave_number-1]
z = wn_data[meridional_mode_type][vertical_wave_numebr-1][2].T[meridional_wave_number-1]

plt.plot(u, lat, label='u')
plt.plot(v, lat, label='v')
plt.plot(z, lat, label='z')
# make a dashed line at 0
plt.axvline(x=0, color='cyan', linestyle='--')

plt.title(f'UVZ profile for {target_hj} : \nzonal k = {zonal_wave_number}, meridional l = {meridional_wave_number}, vertical n = {vertical_wave_numebr}')
plt.ylabel('Latitude')
plt.ylim(-90,90)
plt.legend()
plt.tight_layout()
plt.show()


wn_data.close()