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
          'figure.figsize': (16, 9),
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

my_dir = './SPS.USER.LHCMD4-MD_CRAB_26_200_L3034_Q26_2021_V1/'
coasts = ['coast3_setting1', 'coast3_setting2', 'coast3_setting3']#, 'coast1_setting3', 'coast1_setting3b']
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0, (3, 5, 1, 5, 1, 5))]
my_set = 1
t_corr = 2*3600 # correction for acquisition time
fig, ax = plt.subplots(1)

for i, coast in enumerate(coasts):
    path2files = f'{my_dir}{coast}/'
    files_list = sorted(glob.glob(path2files+'*PM1*')) # we will always use "PM1" (in fact we will set the wirescanners such that PM1 will be the best).
# The feature of the “Best channel” from the wirescanner firmware does not work correctly yet (23Sep2021)

    files2ignore_list = [f'{path2files}2021.10.13.12.15.24.135000_SPS.BWS.41677.V-PM1.parquet', f'{path2files}2021.10.13.12.15.24.135000_SPS.BWS.51637.H-PM1.parquet']

    x_dict, y_dict = {}, {}

    x_dict[f'emit_set{my_set}'], y_dict[f'emit_set{my_set}'] = [], []
    x_dict[f'days_set{my_set}'], y_dict[f'days_set{my_set}'] = [], []
    for filename in files_list:
        if filename not in files2ignore_list:
            print(filename)
            # Load data 
            data = ds.parquet_to_awkward(filename) # type: awkward.highlevel.Array
            #print(data.fields) # print the keys of the awkward array

            acq = data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6 # sec

            if '51637.H' in filename:
                x_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0])
                x_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
            if '41677.V' in filename:
                y_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0])
                y_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
        else:
            print(f'file {filename} ignored')



    # compute the emit grwoth in m/day, only for set 1 and only for set 2
    [mX, bX], covX = np.polyfit(x_dict[f'days_set{my_set}'], x_dict[f'emit_set{my_set}'], deg=1, cov=True)
    [mY, bY], covY = np.polyfit(y_dict[f'days_set{my_set}'], y_dict[f'emit_set{my_set}'], deg=1, cov=True)
    
    errX = np.sqrt(np.diag(covX))[0]
    errY = np.sqrt(np.diag(covY))[0]


    xfmt = md.DateFormatter('%H:%M:%S')

    ax.plot(x_dict[f'days_set{my_set}'], np.array(x_dict[f'emit_set{my_set}'])*1e6, 'o', c='b')
    ax.plot(y_dict[f'days_set{my_set}'], np.array(y_dict[f'emit_set{my_set}'],)*1e6, 'o', c='r')

    ax.plot(x_dict[f'days_set{my_set}'], (np.array(x_dict[f'days_set{my_set}'])*mX+bX)*1e6, c='b', ls=linestyles[i], label=r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX*1e6/24:.2f}'+r'$\pm$'+f'{errX*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    ax.plot(y_dict[f'days_set{my_set}'], (np.array(y_dict[f'days_set{my_set}'])*mY+bY)*1e6, c='r', ls=linestyles[i], label=r'$\mathrm{d\epsilon_y/dt}$'+f'= {mY*1e6/24:.2f}'+r'$\pm$'+f'{errY*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

    plt.ylim(1.0, 10.0)
    
    ax.set_title(f'Coast3, Set {my_set}')
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$\mathrm{\epsilon \ [\mu m]}$')
    ax.legend(loc=2, frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.gca().xaxis.set_major_formatter(xfmt)
    plt.grid(ls='--')
    plt.tight_layout()


plt.savefig(f'emit_grwoth_vs_time_coast3_set{my_set}.png')
