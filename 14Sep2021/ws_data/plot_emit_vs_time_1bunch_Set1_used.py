import os
import datascout as ds
import awkward as ak
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as md
import glob 

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

#path2files = './example_WS_data_19Aug2021/'
path2files = './SPS.USER.LHCMD4-MD_CRAB_26_200_L3034_Q26_2021_V1/'

files_list = sorted(glob.glob(path2files+'*PM2*'))
print(files_list)

# files to be ignored from the plotting and the fit
files2ignore_list = ['2021.08.19.13.49.24.135000_SPS.BWS.41677.V-PM1.parquet']
files2ignore_list = []

entry = 0
subentries = np.arange(0,80) # how many bunches
subsubentry=0

# select bunch
bunch = subentries[0]
# select IN (Set 1) or OUT (Set 2)
my_set = 'Set1' # 'Set2'

emitH_list, emitV_list, daysH_list, daysV_list  = [], [], [], []

for filename in files_list:
    if filename not in files2ignore_list:
        data = ds.parquet_to_awkward(filename) # type: awkward.highlevel.Array
        pd_data = ak.to_pandas(data) # convert awkward arrays to pandas for easier manipilation

        acq = pd_data['cycleStamp'][0][bunch][0]/1e9+ pd_data[f'acq_time_{my_set}'][0][bunch][0]/1e6 # sec
        t_corr = 2*3600


        if '51637.H' in filename:
            emitH_list.append(pd_data[f'emittance_{my_set}'][0][bunch][0])
            daysH_list.append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
        if '41677.V' in filename:
            emitV_list.append(pd_data[f'emittance_{my_set}'][0][bunch][0])
            daysV_list.append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
    else:
        print(f'file {filename} ignored')



# compute the emit grwoth in m/day
[mH, bH], covH = np.polyfit(daysH_list, emitH_list, deg=1, cov=True)
[mV, bV], covV = np.polyfit(daysV_list, emitV_list, deg=1, cov=True)

errH = np.sqrt(np.diag(covH))[0]
errV = np.sqrt(np.diag(covV))[0]

print ("Slope H: " + str(mH))
print ("Intercept : " + str(bH))

print ("Slope V: " + str(mV))
print ("Intercept : " + str(bV))

fig, ax = plt.subplots(1)

xfmt = md.DateFormatter('%H:%M:%S')

ax.plot(daysH_list, np.array(emitH_list)*1e6, 'o', c='b')
ax.plot(daysV_list, np.array(emitV_list)*1e6, 'o', c='r')

ax.plot(daysH_list, (np.array(daysH_list)*mH+bH)*1e6, c='b', label=r'$\mathrm{d\epsilon_H/dt}$'+f'= {mH*1e6/24:.2f}'+r'$\pm$'+f'{errH*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
ax.plot(daysV_list, (np.array(daysV_list)*mV+bV)*1e6, c='r', label=r'$\mathrm{d\epsilon_V/dt}$'+f'= {mV*1e6/24:.2f}'+r'$\pm$'+f'{errV*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

ax.set_title(f'bunch {bunch}, {my_set}')
ax.set_xlabel('Time')
ax.set_ylabel(r'$\mathrm{\epsilon \ [\mu m]}$')
ax.legend()
plt.setp(ax.get_xticklabels(), rotation=45)
plt.gca().xaxis.set_major_formatter(xfmt)
plt.grid(ls='--')
plt.tight_layout()

#plt.savefig(f'emit_vs_time_bunch{bunch}_{my_set}.png', bbox_inches='tight')

plt.show()
plt.close()

