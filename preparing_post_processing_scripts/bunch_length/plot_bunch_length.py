import os
import datascout as ds
import awkward as ak
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as md

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



path2files = '/home/natalia/PhD_projects/CC_MD_2021/abwlm_data/14Sep2021/SPS.USER.LHCMD4/'

files_list = sorted(os.listdir(path2files)[:-1])
print(files_list)

# files to be ignored from the plotting and the fit
files2ignore_list =  []

entry = 0
subentries = np.arange(0,80) # how many bunches
subsubentry=0

# select bunch
bunch = subentries[0]
# select IN (Set 1) or OUT (Set 2)
my_set = 'Set1' # 'Set2'

emit_list, days_list  = [], []

for filename in files_list:
    if filename not in files2ignore_list:
        #data = ds.parquet_to_awkward(path2files+filename) # type: awkward.highlevel.Array
        data = ds.parquet_to_dict(path2files+filename)
        #print(data)
        #quit()
        pd_data = pd.DataFrame.from_dict(data)
        #pd_data = ak.to_pandas(data) # convert awkward arrays to pandas for easier manipilation
        print(pd_data)
        quit()

        acq = pd_data['cycleStamp'][0][bunch][0]/1e9+ pd_data[f'acq_time_{my_set}'][0][bunch][0]/1e6 # sec
        t_corr = 2*3600
        days_list.append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.

        emit_list.append(pd_data[f'bunchLengthMean'][0][bunch][0])
    else:
        print(f'file {filename} ignored')



# compute the emit grwoth in m/day
[m, b], cov = np.polyfit(days_list, emit_list, deg=1, cov=True)
err = np.sqrt(np.diag(cov))[0]
print ("Slope : " + str(m))
print ("Intercept : " + str(b))


fig, ax = plt.subplots(1)

xfmt = md.DateFormatter('%H:%M:%S')

ax.plot(days_list, np.array(emit_list)*1e6, 'o', c='b')

ax.plot(days_list, (np.array(days_list)*m+b)*1e6, c='k', label=r'$\mathrm{d\epsilon/dt}$'+f'= {m*1e6/24:.2f}'+r'$\pm$'+f'{err*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

ax.set_title(f'bunch {bunch}, {my_set}')
ax.set_xlabel('Time')
ax.set_ylabel(r'$\mathrm{\epsilon \ [\mu m]}$')
ax.legend()
plt.setp(ax.get_xticklabels(), rotation=45)
plt.gca().xaxis.set_major_formatter(xfmt)
plt.grid(ls='--')
plt.tight_layout()

plt.savefig(f'emit_vs_time_bunch{bunch}_{my_set}.png', bbox_inches='tight')

plt.show()
plt.close()

