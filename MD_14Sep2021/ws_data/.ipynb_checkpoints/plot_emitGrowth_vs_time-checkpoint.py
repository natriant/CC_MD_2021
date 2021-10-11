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
                x_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0])
                x_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
            if '41677.V' in filename:
                y_dict[f'emit_set{my_set}'].append(data[f'emittance_Set{my_set}'][0][0])
                y_dict[f'days_set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
        else:
            print(f'file {filename} ignored')
            
# compute the emit grwoth in m/day, only for set 1 and only for set 2

for my_set in range(1,3):
    [mX, bX], covX = np.polyfit(x_dict[f'days_set{my_set}'], x_dict[f'emit_set{my_set}'], deg=1, cov=True)
    [mY, bY], covY = np.polyfit(y_dict[f'days_set{my_set}'], y_dict[f'emit_set{my_set}'], deg=1, cov=True)
    
    errX = np.sqrt(np.diag(covX))[0]
    errY = np.sqrt(np.diag(covY))[0]

    print ("Slope X: " + str(mX))
    print ("Intercept : " + str(bX))

    print ("Slope Y: " + str(mY))
    print ("Intercept : " + str(bY))

    
    fig, ax = plt.subplots(1)

    xfmt = md.DateFormatter('%H:%M:%S')

    ax.plot(x_dict[f'days_set{my_set}'], np.array(x_dict[f'emit_set{my_set}'])*1e6, 'o', c='b')
    ax.plot(y_dict[f'days_set{my_set}'], np.array(y_dict[f'emit_set{my_set}'],)*1e6, 'o', c='r')

    ax.plot(x_dict[f'days_set{my_set}'], (np.array(x_dict[f'days_set{my_set}'])*mX+bX)*1e6, c='b', label=r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX*1e6/24:.2f}'+r'$\pm$'+f'{errX*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    ax.plot(y_dict[f'days_set{my_set}'], (np.array(y_dict[f'days_set{my_set}'])*mY+bY)*1e6, c='r', label=r'$\mathrm{d\epsilon_y/dt}$'+f'= {mY*1e6/24:.2f}'+r'$\pm$'+f'{errY*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

    plt.ylim(1.00, 2.6)
    
    ax.set_title(f'Set {my_set}')
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$\mathrm{\epsilon \ [\mu m]}$')
    ax.legend(loc=4)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.gca().xaxis.set_major_formatter(xfmt)
    plt.grid(ls='--')
    plt.tight_layout()

    plt.savefig(f'emit_vs_time_Set{my_set}.png', bbox_inches='tight')
    
## compute average between Set 1 and Set 2
emitX_list_mean = [statistics.mean(k) for k in zip(x_dict[f'emit_set1'], x_dict[f'emit_set2'])]
emitX_list_std =  [statistics.stdev(k) for k in zip(x_dict[f'emit_set1'], x_dict[f'emit_set2'])]
daysX_list_mean = [statistics.mean(k) for k in zip(x_dict[f'days_set1'], x_dict[f'days_set2'])]

emitY_list_mean = [statistics.mean(k) for k in zip(y_dict[f'emit_set1'], y_dict[f'emit_set2'])]
emitY_list_std =  [statistics.stdev(k) for k in zip(y_dict[f'emit_set1'], y_dict[f'emit_set2'])]
daysY_list_mean = [statistics.mean(k) for k in zip(y_dict[f'days_set1'], y_dict[f'days_set2'])]



## Compute the emit grwoth in m/day, averaging Set 1 and Set 2

# Mean Set 1 and Set 2
##### WEIGHTED POLYFIT ######
# documentation: https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
# what does it mean for gaussian uncertainties --> sigma is the error of the value
weightsX = 1/np.array(emitX_list_std)
weightsY = 1/np.array(emitY_list_std)
#Alternative (?) weights = 1/np.array(emit_list_std[1:])**2
[mX_mean, bX_mean], covX = np.polyfit(daysX_list_mean, emitX_list_mean, deg=1, w=weightsX, cov=True)
errX_mean = np.sqrt(np.diag(covX))[0]
[mY_mean, bY_mean], covY = np.polyfit(daysY_list_mean, emitY_list_mean, deg=1, w=weightsY, cov=True)
errY_mean = np.sqrt(np.diag(covY))[0]


print ("Slope x-plane: " + str(mY_mean))
print ("Intercept x-plane: " + str(bX_mean))
print ("Slope y-plane: " + str(mY_mean))
print ("Intercept y-plane: " + str(bX_mean))

fig, ax = plt.subplots(1)

xfmt = md.DateFormatter('%H:%M:%S')

ax.errorbar(daysX_list_mean, np.array(emitX_list_mean)*1e6, yerr=np.array(emitX_list_std)*1e6, marker='o', ls='',c='b', capsize=5, zorder=50, label=r'$\mathrm{\langle Set 1, Set 2 \rangle}$')
ax.plot(daysX_list_mean, (np.array(daysX_list_mean)*mX_mean+bX_mean)*1e6, c='b', label=r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX_mean*1e6/24:.2f}'+r'$\pm$'+f'{errX_mean*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

ax.errorbar(daysY_list_mean, np.array(emitY_list_mean)*1e6, yerr=np.array(emitY_list_std)*1e6, marker='o', ls='',c='r', capsize=5, zorder=50, label=r'$\mathrm{\langle Set 1, Set 2 \rangle}$')
ax.plot(daysY_list_mean, (np.array(daysY_list_mean)*mY_mean+bY_mean)*1e6, c='r', label=r'$\mathrm{d\epsilon_y/dt}$'+f'= {mY_mean*1e6/24:.2f}'+r'$\pm$'+f'{errY_mean*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')

plt.ylim(1.00, 2.6)

ax.set_xlabel('Time')
ax.set_ylabel(r'$\mathrm{\epsilon \ [\mu m]}$')
ax.legend(loc=4)
plt.setp(ax.get_xticklabels(), rotation=45)
plt.gca().xaxis.set_major_formatter(xfmt)
plt.grid(ls='--')
plt.tight_layout()


plt.savefig(f'emit_vs_time_averageSet1andSet2.png', bbox_inches='tight')
