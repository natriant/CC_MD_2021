import pickle
import numpy as np
from datetime import datetime
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
plt.ion()

mr_overview = pickle.load(open('mr_overview.pkl', 'rb'))
print(mr_overview['LHCMD4'].keys())

bunch_length = mr_overview['LHCMD4']['sigma_t']
timestamps = mr_overview['LHCMD4']['timestamp_float']

t_corr = 2*3600
acq = timestamps
my_title = str(datetime.fromtimestamp(timestamps[0]))[:10] # Y-M-D

my_time = md.epoch2num(acq+t_corr)[120:287] # Convert UNIX time to days since Matplotlib epoch.
bunch_length_cut = bunch_length[120:287] # [120:290] Acquisitions that correspond to the time of the emittance growth measurements. [120:287] chosen to analyse exactly 1h of measurements.
# Index : 287 to be exactly 1h of measurements (real time up to 290)

fig, ax = plt.subplots(1)

xfmt = md.DateFormatter('%H:%M:%S')
ax.plot(my_time, bunch_length_cut*4, 'o', c='lightgreen')
ax.set_ylim(1.9, 2.5)
ax.set_xlabel('Time')
ax.set_ylabel('$\mathrm{4 \sigma_t \ [ns]}$')
ax.set_title(f'Mountain Range: {my_title}')

plt.setp(ax.get_xticklabels(), rotation=45)
plt.gca().xaxis.set_major_formatter(xfmt)
plt.grid(ls='--')
plt.tight_layout()

savefig = False
if savefig:
    plt.savefig(f'./figures/MR_{my_title}.png', bbox_inches='tight')
                                   
                                   
# Linear fit on the bunch length evolution
[m, b], cov = np.polyfit(my_time, bunch_length_cut*4, deg=1, cov=True)
err = np.sqrt(np.diag(cov))[0]
print(m, b, cov)
                                   
# Compute bunch lenth increase in percentage. 
# As there is a flactuation in the bunch length values, use the average of the 5 first and 5 lasts measurements.
bunch_length_init, bunch_length_fin = np.mean(bunch_length_cut[:5]), np.mean(bunch_length_cut[-5:])  
D_sigma = bunch_length_fin-bunch_length_init
perc = (D_sigma/bunch_length_init)*100

# Plot bunch length evolution with linear fit                                   
fig, ax = plt.subplots(1)

xfmt = md.DateFormatter('%H:%M:%S')
ax.plot(my_time, bunch_length_cut*4, 'o', c='lightgreen')
ax.plot(my_time, b+m*my_time, c='k', label=r'$\mathrm{4 d\sigma_t/dt}$'+f'= {m/24:.2f}'+r'$\mathrm{[ns/h]}$'+r'$\equiv $'+f'{perc:.2f}'+r'$\%$'+'/h')

ax.set_ylim(1.9, 2.5)
ax.set_xlabel('Time')
ax.set_ylabel('$\mathrm{4 \sigma_t \ [ns]}$')
ax.set_title(f'Mountain Range: {my_title}')
ax.legend(loc=3)

plt.setp(ax.get_xticklabels(), rotation=45)
plt.gca().xaxis.set_major_formatter(xfmt)
plt.grid(ls='--')
plt.tight_layout()

savefig = True
if savefig:
    plt.savefig(f'./figures/MR_{my_title}_fit.png', bbox_inches='tight')
    

