import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import netCDF4 as nc

# =======================
PLOT_VSF = False
PLOT_HSF = False
PRINT_DATA = False


# Read in the data
filepath_freq = './hot_jupiter/MODES/Data/Output for HD 209458 b/freq_0707.dat'
filename_freq = os.path.basename(filepath_freq)
filepath_wn = '/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for HD 209458 b/wn_folder/hough.wn00000.nc'
wn_0 = nc.Dataset(filepath_wn)
filepath_vsf = '/Users/changyichieh/Documents/Hot Jupiter/Hot_Jupiter_Github/hot_jupiter/MODES/Data/Output for Earth/vsf.data_Earth.nc'
vsf = nc.Dataset(filepath_vsf)



def plot_vsf_modes():
    # plot the vertical structure function for the first 10 modes
    plt.figure(figsize=(10, 8))
    plt.title('Vertical Structure Function for the first 5 modes : HD 209458 b')
    plt.xlabel('Vertical structure function')
    plt.ylabel('Sigma')
    plt.gca().invert_yaxis()
    equi_height = vsf['evht'][:]
    # print(vsf['vgrid_weight'][:])
    plt.plot(vsf['vsf'][0].data, vsf['vgrid'][:].data, label=f'Mode 0, Equivalent Height = {equi_height[0]:.2f}')
    plt.plot(vsf['vsf'][1].data, vsf['vgrid'][:].data, label=f'Mode 1, Equivalent Height = {equi_height[1]:.2f}')
    plt.plot(vsf['vsf'][2].data, vsf['vgrid'][:].data, label=f'Mode 2, Equivalent Height = {equi_height[2]:.2f}')
    plt.plot(vsf['vsf'][3].data, vsf['vgrid'][:].data, label=f'Mode 3, Equivalent Height = {equi_height[3]:.2f}')
    plt.plot(vsf['vsf'][4].data, vsf['vgrid'][:].data, label=f'Mode 4, Equivalent Height = {equi_height[4]:.2f}')
    # Draw a vertical line at 0
    plt.axvline(x=0, color='k', linestyle='--')
    plt.legend()
    plt.xlim(-0.4,0.4)
    plt.yscale('log')
    plt.show()

def get_evht(data_path):
    data = nc.Dataset(data_path)
    evht = data['evht'][:]
    return evht



if PLOT_VSF:
# plot the first 5 vertical modes (VSF solutions)
    plot_vsf_modes()

def process_freq(dataset_path):

    # count the number of lines in the file (each of the data is stored in 23 lines)
    num_data = sum(1 for line in open(dataset_path))/23

    # Create a empty dataframe and assign the header
    df = pd.DataFrame(columns=['Vertical_mode_number', 'eps', 'Equivalent Height', 'zonal_wave_number', 'Eigenfrequency of eastward gravity mode', 'Eigenfrequency of westward gravity mode', 'Eigenfrequency of rotaitonal mode'])

    # Read in the data and store them in a list for every 23 lines
    with open(dataset_path,'r') as f:
        data = f.readlines()

        for i in range(int(num_data)):
            # for "Vertical_mode_number     eps     Equivalent Height zonal_wave_number"
            float_list = [float(value) for value in data[i*23+1].split()]
            # for "Eigenfrequency of eastward gravity mode"
            data_east_gravity = data[i*23+3:i*23+9]
            Eastward_gravity_list = []

            for line in data_east_gravity:
                values = line.split()
                float_values = [float(value) for value in values]
                Eastward_gravity_list.extend(float_values)

          # for "Eigenfrequency of westward gravity mode"
            data_west_gravity = data[i*23+10:i*23+16]
            Westward_gravity_list = []

            for line in data_west_gravity:
                values = line.split()
                float_values = [float(value) for value in values]
                Westward_gravity_list.extend(float_values)


            # for "Eigenfrequency of rotaitonal mode"
            data_rotational = data[i*23+17:i*23+23]
            Rotational_list = []

            for line in data_rotational:
                values = line.split()
                float_values = [float(value) for value in values]
                Rotational_list.extend(float_values)

            # store the above data to the dataframe   
            df = df.append({'Vertical_mode_number': float_list[0], 
                    'eps': float_list[1], 
                    'Equivalent Height': float_list[2], 
                    'zonal_wave_number': float_list[3], 
                    'Eigenfrequency of eastward gravity mode': Eastward_gravity_list, 
                    'Eigenfrequency of westward gravity mode': Westward_gravity_list, 
                    'Eigenfrequency of rotaitonal mode': Rotational_list}, ignore_index=True)
    return df
    # print(type(data[1]))

df = process_freq(filepath_freq)

# clear the output terminal
os.system('clear')

# print the headers of the dataframe
headers = df.columns.tolist()
print('-------------------headers-------------------')
for header in headers:
    print(header)
print('---------------------------------------------')

# print the data
if PRINT_DATA:
    print(df['Eigenfrequency of eastward gravity mode'])
    print('----------')
    print(df['Eigenfrequency of westward gravity mode'])
    print('----------')
    print(df['Eigenfrequency of rotaitonal mode'])

#===================================================================================================
# enter the desired wavenumber
# !! Notice that the wavenumber of zonal wave starts from 0, and the wavenumber of meridional/verical starts from 1
zonal_wave_number = 1           # 0 means the "first" zonal mode but the index is 0 ; 41 in total (0~40)
vertical_wave_number = 5       # 1 means the first vertical mode ; 0 means nothing ; 20 in total 
meridional_wave_type = 'EIG'    # 'WIG' means westward inertial gravity mode ; 'EIG' means eastward inertial gravity mode ; 'ROT' means rotational mode each has 30 modes
meridional_wave_number = 1      # 1 means the first meridional mode ; 0 means nothing ; 30 in total
#===================================================================================================

def find_eigenfrequncy(zonal_wave_number, vertical_wave_number, meridional_wave_type, meridional_wave_number):
    eigenfrequency = -100000
    # find the index of the desired wavenumber
    index = df.index[(df['zonal_wave_number'] == zonal_wave_number) & (df['Vertical_mode_number'] == vertical_wave_number)].tolist()
    result = df.iloc[index]

    # find the eigenfrequncy with the desired meridional wave number
    if meridional_wave_type == 'WIG':
        eigenfrequency = result['Eigenfrequency of westward gravity mode'].values[0][meridional_wave_number-1]
    elif meridional_wave_type == 'EIG':
        eigenfrequency = result['Eigenfrequency of eastward gravity mode'].values[0][meridional_wave_number-1]
    elif meridional_wave_type == 'ROT':
        eigenfrequency = result['Eigenfrequency of rotaitonal mode'].values[0][meridional_wave_number-1]
    return eigenfrequency

# # check the result for (0,0,1) to (0,20,1) (zonal, vertical, meridional) [in the note, (0,20,1) means (0,19,0) here]
# for i in range(20):
#     eigenfrequency = find_eigenfrequncy(zonal_wave_number, i+1, meridional_wave_type, meridional_wave_number)
#     print(f'The Eigenfrequency of zonal wave number {zonal_wave_number}, vertical wave number {i}, meridional type {meridional_wave_type} wave number {meridional_wave_number} : {eigenfrequency}')

# eigenfrequency = find_eigenfrequncy(zonal_wave_number, vertical_wave_number, meridional_wave_type, meridional_wave_number)
# # print the result
# print(f'The Eigenfrequency of zonal wave number {zonal_wave_number}, vertical wave number {vertical_wave_number}, meridional type {meridional_wave_type} wave number {meridional_wave_number} : {eigenfrequency}')



def plot_hsf(zonal_wave_number,meridional_wave_number,meridional_wave_type):
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
        plt.plot(sqrt_epsilon, frequency_list, label=rf'{meridional_wave_type} $\ell=${i}')
    plt.xlabel(r'$\sqrt{\epsilon}$')
    plt.ylabel('eigenfrequency')
    # if meridional_wave_type != 'EIG':
    #     plt.gca().invert_yaxis()
    plt.gca().invert_yaxis()
    plt.title(f'{meridional_wave_type} mode, zonal wave number k={zonal_wave_number}')
    plt.yscale('symlog', linthresh=0.0001)
    plt.xscale('log')
    plt.xlim(0.1,10)
    # plt.ylim(-0.001,-1)
    plt.grid(True)
    plt.legend()
    plt.show()





# exit()
if PLOT_HSF: 
    plot_hsf(zonal_wave_number,meridional_wave_number,meridional_wave_type)

# close the file
wn_0.close()
vsf.close()