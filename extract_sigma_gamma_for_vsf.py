import pandas as pd
import numpy as np
import os

# Read in the data
filepath = 'hot_jupiter/Parmentier 2/output/PTprofile(20bar)(Hd209458b)(1107).csv'
filename = os.path.basename(filepath)
df = pd.read_csv(filepath)
df = df[:-1]

print('----------------------File Name----------------------')
print(filename)
print('-----------------------------------------------------')

# Extract the data
sigma = df['sigma'].values
gamma = df['GAMMA_0'].astype(float).values

# reverse the order of the data
sigma = sigma[::-1]
gamma = gamma[::-1]

# Write the data into a new .txt file seperately
with open('hot_jupiter/Parmentier 2/output/sigma.txt', 'w') as f:
    for item in sigma:
        f.write("%s\n" % item)

with open('hot_jupiter/Parmentier 2/output/gamma.txt', 'w') as f:
    for item in gamma:
        f.write("%s\n" % item)


