import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import netCDF4 as nc
from analyzing_hj_mode import process_freq, get_evht
import sys


#===================================================================================================
# enter the desired wavenumber
# !! Notice that the wavenumber of zonal wave starts from 0, and the wavenumber of meridional/verical starts from 1
zonal_wave_number = 1           # 0 means the "first" zonal mode but the index is 0 ; 41 in total (0~40)
# vertical_wave_number = 20       # 1 means the first vertical mode ; 0 means nothing ; 20 in total 
vertical_wave_number = int(sys.argv[1])
meridional_wave_type = 'EIG'    # 'WIG' means westward inertial gravity mode ; 'EIG' means eastward inertial gravity mode ; 'ROT' means rotational mode each has 30 modes
meridional_wave_number = 3     # 1 means the first meridional mode ; 0 means nothing ; 30 in total
PLOT_MERI = False
PLOT_EVHT = False
PLOT_HOUGH = True
#===================================================================================================
# Assuming you have meridional_wave_number as the maximum value for normalization
normalize = plt.Normalize(vmin=0, vmax=meridional_wave_number)


# for WASP19b
filepath_freq_WASP = 'MODES/Data/Output for WASP 19 b/freq_0718.dat'
filepath_vsf_WASP = '/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for WASP 19 b/vsf.data.nc'
df_WASP = process_freq(filepath_freq_WASP)
df_WASP_zonal = df_WASP[df_WASP['zonal_wave_number'] == zonal_wave_number]
df_WASP_zonal = df_WASP_zonal[df_WASP_zonal['Vertical_mode_number'] <= vertical_wave_number]

sqrt_WASP_epsilon = df_WASP_zonal['eps']**0.5
if meridional_wave_type == 'WIG':
    df_WASP_hsf = df_WASP_zonal['Eigenfrequency of westward gravity mode']
elif meridional_wave_type == 'EIG':
    df_WASP_hsf = df_WASP_zonal['Eigenfrequency of eastward gravity mode']
elif meridional_wave_type == 'ROT':
    df_WASP_hsf = df_WASP_zonal['Eigenfrequency of rotaitonal mode']
df_WASP_hsf = df_WASP_hsf.values.tolist()

# Replace 'cool' with any other valid colormap name from Matplotlib
cmap_WASP = plt.get_cmap('cool')

if PLOT_HOUGH:
    df_WASP_Wn = nc.Dataset(f"MODES/Data/Output for WASP 19 b/hough.nc/hough.wn{zonal_wave_number:05d}.nc")
    lat = df_WASP_Wn['lat'][:]
    if meridional_wave_type == 'ROT': meridional_wave_type = 'BAL' # change the name due to the data stucture
    u = df_WASP_Wn[meridional_wave_type][vertical_wave_number-1][0].T[meridional_wave_number-1]
    v = df_WASP_Wn[meridional_wave_type][vertical_wave_number-1][1].T[meridional_wave_number-1]
    z = df_WASP_Wn[meridional_wave_type][vertical_wave_number-1][2].T[meridional_wave_number-1]
    plt.figure(figsize=(6,9))
    
    plt.subplot(1,2,1)
    plt.axvline(x=0, color='k', linestyle='--',alpha=0.5)
    plt.plot(u, lat, label='u')
    plt.plot(v, lat, label='v')
    plt.title('WASP19b')
    plt.ylabel('Latitude')
    plt.xlim(-1,1)
    plt.ylim(-90,90)
    plt.legend()




if PLOT_MERI:
    for i in range(meridional_wave_number):
        frequency_WASP_list = []
        
        for j in range(len(sqrt_WASP_epsilon)):
            # frequency_list.append(abs(df_hsf[j][i])) # abs() is to make the negative value positive
            frequency_WASP_list.append(df_WASP_hsf[j][i]) 
        # plt.scatter(sqrt_WASP_epsilon, frequency_WASP_list,label=rf'{meridional_wave_type} $\ell=${i} for WASP19b')
        color = cmap_WASP(normalize(i))
        plt.plot(sqrt_WASP_epsilon, frequency_WASP_list, marker='s',markersize=2,label=rf'{meridional_wave_type} $\ell=${i} for WASP19b',color=color)

if PLOT_EVHT:
    evht_WASP = get_evht(filepath_vsf_WASP)
    number_WASP = np.ones(len(evht_WASP))*1
    # plt.scatter(number_WASP,evht_WASP,  marker='.',color='blue',label='WASP19b',s=5)
    plt.plot(number_WASP,evht_WASP,  marker='.',color='blue',linestyle='dashed',markersize=5,label='WASP19b')

# for Hd209458b
filepath_freq_HD = 'MODES/Data/Output for HD 209458 b/freq_0707.dat'
filepath_vsf_HD = '/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for HD 209458 b/vsf.data.nc'
df_HD = process_freq(filepath_freq_HD)
df_HD_zonal = df_HD[df_HD['zonal_wave_number'] == zonal_wave_number]
df_HD_zonal = df_HD_zonal[df_HD_zonal['Vertical_mode_number'] <= vertical_wave_number]
sqrt_HD_epsilon = df_HD_zonal['eps']**0.5
if meridional_wave_type == 'WIG':
    df_HD_hsf = df_HD_zonal['Eigenfrequency of westward gravity mode']
elif meridional_wave_type == 'EIG':
    df_HD_hsf = df_HD_zonal['Eigenfrequency of eastward gravity mode']
elif meridional_wave_type == 'ROT':
    df_HD_hsf = df_HD_zonal['Eigenfrequency of rotaitonal mode']
df_HD_hsf = df_HD_hsf.values.tolist()
cmap_HD = plt.get_cmap('hot')

if PLOT_HOUGH:
    df_HD_Wn = nc.Dataset(f"/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for HD 209458 b/wn_folder/hough.wn{zonal_wave_number:05d}.nc")
    lat = df_HD_Wn['lat'][:]
    if meridional_wave_type == 'ROT': meridional_wave_type = 'BAL' # change the name due to the data stucture
    u = df_HD_Wn[meridional_wave_type][vertical_wave_number-1][0].T[meridional_wave_number-1]
    v = df_HD_Wn[meridional_wave_type][vertical_wave_number-1][1].T[meridional_wave_number-1]
    z = df_HD_Wn[meridional_wave_type][vertical_wave_number-1][2].T[meridional_wave_number-1]
    plt.subplot(1,2,2)
    plt.axvline(x=0, color='k', linestyle='--',alpha=0.5)
    plt.plot(u, lat, label='u')
    plt.plot(v, lat, label='v')
    plt.title('Hd209458b')
    plt.ylabel('Latitude')
    plt.xlim(-1,1)
    plt.ylim(-90,90)
    plt.legend()
    plt.suptitle(f'zonal k = {zonal_wave_number}, vertical n = {vertical_wave_number}, meridional {meridional_wave_type} m = {meridional_wave_number}')
    plt.savefig(dpi=300,fname=f'zonal{zonal_wave_number:02d}_vertical{vertical_wave_number:02d}_meridional{meridional_wave_type}_{meridional_wave_number:02d}.png')
    plt.tight_layout()
    


if PLOT_MERI:
    for i in range(meridional_wave_number):
        frequency_HD_list = []
        for j in range(len(sqrt_HD_epsilon)):
            # frequency_list.append(abs(df_hsf[j][i])) # abs() is to make the negative value positive
            frequency_HD_list.append(df_HD_hsf[j][i]) 
        color = cmap_HD(normalize(i))
        plt.plot(sqrt_HD_epsilon, frequency_HD_list,linestyle='dashed', marker='.',markersize=5,label=rf'{meridional_wave_type} $\ell=${i} for Hd209458b',color=color,alpha=0.5)

if PLOT_EVHT:
    evht_HD = get_evht(filepath_vsf_HD)
    number_HD = np.ones(len(evht_WASP))*2
    # plt.scatter( number_HD,evht_HD, marker='.',color='green',label='Hd209458b',s=5)
    plt.plot( number_HD,evht_HD, marker='.',color='green',linestyle='dashed',markersize=5,label='Hd209458b')



if PLOT_EVHT:
    evht_Earth = get_evht('/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for Earth/vsf.data_Earth.nc')
    number_Earth = np.ones(len(evht_Earth))*3
    plt.plot( number_Earth,evht_Earth, marker='.',color='orange',linestyle='dashed',markersize=5,label='Earth')
    # plt.yscale('log')
    print(evht_HD[0])
    print(evht_WASP[0])
    plt.ylabel('Equivalent Height')
    
    plt.legend()
    plt.xlim(0,4)
    plt.show()


if PLOT_MERI:
    plt.xlabel(r'$\sqrt{\epsilon}$')
    plt.ylabel('eigenfrequency')
    # if meridional_wave_type != 'EIG':
    #     plt.gca().invert_yaxis()

    plt.title(f'{meridional_wave_type} mode, zonal wave number k={zonal_wave_number}')
    plt.yscale('symlog', linthresh=0.0001)
    plt.xscale('log')
    # plt.xlim(0.1,10)
    if meridional_wave_type == 'WIG':
        plt.ylim(-10,-0.1)
        plt.gca().invert_yaxis()
    elif meridional_wave_type == 'EIG':
        plt.ylim(0.01,10)
    elif meridional_wave_type == 'ROT':
        plt.ylim(-0.001,-1)
        # plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.show()


exit()




def plot_hsf_v2(data,zonal_wave_number,meridional_wave_number,meridional_wave_type):
    df = data
    # create a new datafram containing the desired zonal wave number
    df_zonal = df[df['zonal_wave_number'] == zonal_wave_number]
    print(df_zonal)
    sqrt_epsilon = df_zonal['eps']**0.5
    if meridional_wave_type == 'WIG':
        df_hsf = df_zonal['Eigenfrequency of westward gravity mode']
    elif meridional_wave_type == 'EIG':
        df_hsf = df_zonal['Eigenfrequency of eastward gravity mode']
    elif meridional_wave_type == 'ROT':
        df_hsf = df_zonal['Eigenfrequency of rotaitonal mode']
    df_hsf = df_hsf.values.tolist()
    # df_hsf = pd.DataFrame(df_hsf.values.tolist()).T.values.tolist() # transpose the data
    # df_hsf = pd.DataFrame(df_hsf).abs()
    for i in range(meridional_wave_number):
        frequency_list = []
        for j in range(len(sqrt_epsilon)):
            # frequency_list.append(abs(df_hsf[j][i])) # abs() is to make the negative value positive
            frequency_list.append(df_hsf[j][i]) 

        if data == data_WASP:
            plt.plot(sqrt_epsilon, frequency_list,'-',label=rf'{meridional_wave_type} $\ell=${i} from WASP19b')
        elif data == data_Hd:
            plt.plot(sqrt_epsilon, frequency_list,label=rf'{meridional_wave_type} $\ell=${i} from Hd209468')
        # plt.plot(sqrt_epsilon, frequency_list, label=rf'{meridional_wave_type} $\ell=${i}')
    plt.xlabel(r'$\sqrt{\epsilon}$')
    plt.ylabel('eigenfrequency')
    # if meridional_wave_type != 'EIG':
    #     plt.gca().invert_yaxis()
    plt.gca().invert_yaxis()
    plt.title(f'{meridional_wave_type} mode, zonal wave number k={zonal_wave_number}')
    plt.yscale('symlog', linthresh=0.0001)
    plt.xscale('log')
    # plt.xlim(0.1,10)
    plt.ylim(-0.001,-1)
    plt.grid(True)
    plt.legend()
    # plt.show()

data_WASP = process_freq('MODES/Data/Output for WASP 19 b/freq_0718.dat')
data_Hd = process_freq('MODES/Data/Output for HD 209458 b/freq_0707.dat')

plot_hsf_v2(data_WASP,1,3,'ROT')
plot_hsf_v2(data_Hd,1,3,'ROT')
plt.show()