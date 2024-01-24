import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import check_wn_file
test = nc.Dataset("hot_jupiter/MODES/Data/Output for Earth/vsf.data_Earth.nc")
print("-----------Equivalent Height (M)-----------")
print(test['evht'][:])
print("-----------sqrt(gD)-----------")
print((9.8*test['evht'][:])**0.5)
evht_0 = test['evht'][0]
scale  = (evht_0*9.8)**0.5
test.close()
# exit()


filepath = 'hot_jupiter/MODES/Data/erainterim_moda_L60_N24_200001010000000.nc'
filepath_inverse = 'hot_jupiter/MODES/Data/Output for Earth/Erai_allk_Kelvin_inverse.nc'
filepath_inverse_no_filter_constant_field = 'hot_jupiter/MODES/Data/Output for Earth/Erai_allk_Kelvin_inverse_no_filter_constant_field.nc'
data = nc.Dataset(filepath)
# print(data)
# print(data['lon'][:])
# print(data['lat'][:])
data_inverse = nc.Dataset(filepath_inverse)
data_inverse_no_filter_constant_field = nc.Dataset(filepath_inverse_no_filter_constant_field)

#choose the data from : lon = 0 -> the first 'lon' is 0 ; level = 0 

# print(data['lat'])
print(data_inverse)
# for i in range(len(data_inverse['u'][0][-1].T)):
# for i in 20 to 31

# for i in range(0,len(data_inverse['u'][0][13].T)):
# print(data_inverse['lev'])
    # for j in range(len(data_inverse['u'][0][:])):
lon_index = 48
lat_inverse = data_inverse['lat'][:]
u_inverse = data_inverse['u'][0][0].T[lon_index]*2.8169299359176234
v_inverse = data_inverse['v'][0][0].T[lon_index]*2.8169299359176234
# z_inverse = data_inverse['Z'][0][0].T[48]
print(v_inverse[0])
# print(z_inverse)
# z_inverse /= evht_0
print(-0.03875792771854539/-0.013758925)
# print(z_inverse)
# plt.figure(figsize=(6, 6))
plt.plot(u_inverse, lat_inverse,label='u from input')
plt.plot(v_inverse, lat_inverse,label='v from input',alpha=0.5) 
# plt.plot(z_inverse, lat_inverse,label='z w bk')
plt.title(f'UVZ profile for Earth : \n longtitude = {data_inverse["lon"][lon_index]:.3f}')
# plt.xlim(-0.1,0.8)
plt.legend()
plt.savefig(f"/Users/changyichieh/Desktop/0725temp/result_{47}.png",dpi=300)
print(f"Done {47}")
plt.show()

# lat = data['lat'][:]
# #var129:z;var131:u;var132:v
# u = data['var131'][0][0].T[0]
# v = data['var132'][0][0].T[0]
# z = data['var129'][0][0].T[0]
# # print(u)
# plt.figure(figsize=(6, 10))
# plt.plot(u, lat,label='u from input')
# plt.plot(v, lat,label='v from input')
# # plt.plot(z, lat,label='z from input')
# plt.legend()
# plt.show()

data_inverse_no_filter_constant_field.close()
data_inverse.close()
data.close()