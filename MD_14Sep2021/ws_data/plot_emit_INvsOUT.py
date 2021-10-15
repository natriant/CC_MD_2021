import os
import datascout as ds
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import glob 
import statistics

# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (8, 7),
          'axes.labelsize': 25,
          'axes.titlesize': 21,
          'xtick.labelsize': 23,
          'ytick.labelsize': 23,
          'image.cmap': 'jet',
          'lines.linewidth': 2,
          'lines.markersize': 10,
          'font.family': 'sans-serif'}


plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)
plt.ion()

plt.close('all')

path2files = './SPS.USER.LHCMD4-MD_CRAB_26_200_L3034_Q26_2021_V1/emit_growth/'
files_list = sorted(glob.glob(path2files+'*PM1*')) # we will always use "PM1" (in fact we will set the wirescanners such that PM1 will be the best).
# The feature of the “Best channel” from the wirescanner firmware does not work correctly yet (23Sep2021)

files2ignore_list = []

x_dict, y_dict = {}, {}

t_corr = 2*3600 # correction for acquisition time


for my_set in range(1,3):
    x_dict[f'emit_set{my_set}'], y_dict[f'emit_set{my_set}'] = [], []
    x_dict[f'days_set{my_set}'], y_dict[f'days_set{my_set}'] = [], []
    for filename in files_list:
        if filename not in files2ignore_list:
            # Load data 
            data = ds.parquet_to_awkward(filename) # type: awkward.highlevel.Array
            #print(data.fields) # print the keys of the awkward array

            acq = data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6 # sec

            if '51637.H' in filename:
                x_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0]*1e6)
                x_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
            if '41677.V' in filename:
                y_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0]*1e6)
                y_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
        else:
            print(f'file {filename} ignored')
            

# '51637.H'
fig, ax = plt.subplots()
ax.plot(x_dict[f'emit_set1'], x_dict[f'emit_set2'], '-o', c='b')
ax.set_ylabel(r'$\mathrm{IN \ [\mu m/h]}$') 
ax.set_xlabel(r'$\mathrm{OUT \ [\mu m/h]}$') 
ax.set_title('BWS.51637.H')
plt.grid(ls='--')
plt.tight_layout()

plt.savefig(f'./figures/BWS.51637.H_IN_vs_OUT.png', bbox_inches='tight')


# '41677.V'
fig, ax = plt.subplots()
ax.plot(y_dict[f'emit_set1'], y_dict[f'emit_set2'], '-o', c='r')
ax.set_ylabel(r'$\mathrm{IN \ [\mu m/h]}$') 
ax.set_xlabel(r'$\mathrm{OUT \ [\mu m/h]}$') 
ax.set_title('BWS.41677.V')
plt.grid(ls='--')
plt.tight_layout()

plt.savefig(f'./figures/BWS.41677.V_IN_vs_OUT.png', bbox_inches='tight')