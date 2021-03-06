'''
Compare the emittance values as computed from the profiles (myFit) vs the values from the parquet files.
'''
import os
import datascout as ds
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import glob 
import datetime
import statistics
from utils import *

# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (16, 7),
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
my_variables = ['emittance_Set', 'profiles_Set', 'positions_Set', 'betagamma_Set', 'days_Set']
for my_set in range(1,3):
    for variable in my_variables:
        x_dict[f'{variable}{my_set}'], y_dict[f'{variable}{my_set}'] = [], []
    for filename in files_list:
        if filename not in files2ignore_list:
            # Load data 
            data = ds.parquet_to_awkward(filename) # type: awkward.highlevel.Array
            #print(data.fields) # print the keys of the awkward array
            acq = data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6 # sec

            if '51637.H' in filename:
                x_dict[f'betagamma_Set{my_set}'].append(data[f'betagamma_Set{my_set}'])
                x_dict[f'days_Set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
                #x_dict[f'timestamp_Set{my_set}'].append(data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6)  # Convert UNIX time to days since Matplotlib epoch.
                for variable in my_variables[:-2]:
                    x_dict[f'{variable}{my_set}'].append(data[f'{variable}{my_set}'][0][0])
                   
            if '41677.V' in filename:
                y_dict[f'betagamma_Set{my_set}'].append(data[f'betagamma_Set{my_set}'])
                y_dict[f'days_Set{my_set}'].append(md.epoch2num(acq+t_corr))  # Convert UNIX time to days since Matplotlib epoch.
                #y_dict[f'timestamp_Set{my_set}'].append(data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6)
                for variable in my_variables[:-2]:
                    y_dict[f'{variable}{my_set}'].append(data[f'{variable}{my_set}'][0][0])
        else:
            print(f'file {filename} ignored')

# Optic parameters for Q26 (Hannes)
beta = {'SPS.BWS.51637.H': 80.60665636, 
        'SPS.BWS.41677.V': 61.23085007}

Dx = {'SPS.BWS.51637.H':  -0.1801965987}

dpp = 0.0003865047556 # for Vrf = 3.79 MV and bunch length of 2.2 ns at 200 GeV.


# Perform fit again and add it to the dictionary. Call it my fit or new fit. 
# betagamma variable should be the same for every case. However, here it is loaded from the data every time for adaptability.
for my_set in range(1,3):
    # Intialize lists in the existing directory for the values from the new fit.
    y_dict[f'emittance_Set{my_set}_myFit'] = []
    x_dict[f'emittance_Set{my_set}_myFit'] = []
    x_dict[f'emittance_Set{my_set}_myFit_with_Dx'] = []
    
    for i in range(len(y_dict[f'positions_Set{my_set}'])): # iterate over the acquisitions for each setting, y-plane
        pop_t_y = fitGauss(y_dict[f'positions_Set{my_set}'][i], y_dict[f'profiles_Set{my_set}'][i])
        ey = getEmittance(pop_t_y[2], beta['SPS.BWS.41677.V'], y_dict[f'betagamma_Set{my_set}'][i]) # returns the emittance in [um]
        y_dict[f'emittance_Set{my_set}_myFit'].append(ey[0]*1e-6)
    for i in range(len(x_dict[f'positions_Set{my_set}'])): # iterate over the acquisitions for each setting, x-plane
        pop_t_x = fitGauss(x_dict[f'positions_Set{my_set}'][i], x_dict[f'profiles_Set{my_set}'][i])
        ex = getEmittance(pop_t_x[2], beta['SPS.BWS.51637.H'], x_dict[f'betagamma_Set{my_set}'][i])
        ex_Dx = getEmittance_with_Dx(pop_t_x[2],  beta['SPS.BWS.51637.H'], x_dict['betagamma_Set1'][0], Dx['SPS.BWS.51637.H'], dpp) # returns the emittance in [um]
        x_dict[f'emittance_Set{my_set}_myFit'].append(ex[0]*1e-6) 
        x_dict[f'emittance_Set{my_set}_myFit_with_Dx'].append(ex_Dx[0]*1e-6)


# compute the emit grwoth in m/day, only for set 1 and only for set 2 and plot new vs old
for my_set in range(1,3):
    [mX, bX], covX = np.polyfit(x_dict[f'days_Set{my_set}'], x_dict[f'emittance_Set{my_set}'], deg=1, cov=True)
    [mY, bY], covY = np.polyfit(y_dict[f'days_Set{my_set}'], y_dict[f'emittance_Set{my_set}'], deg=1, cov=True)
    # error on the slop of the fit
    errX = np.sqrt(np.diag(covX))[0]
    errY = np.sqrt(np.diag(covY))[0]

    # myFit
    [mX_myFit, bX_myFit], covX_myFit = np.polyfit(x_dict[f'days_Set{my_set}'], x_dict[f'emittance_Set{my_set}_myFit'], deg=1, cov=True)
    print( x_dict[f'emittance_Set{my_set}_myFit'])
    [mX_myFit_Dx, bX_myFit_Dx], covX_myFit_Dx = np.polyfit(x_dict[f'days_Set{my_set}'], x_dict[f'emittance_Set{my_set}_myFit_with_Dx'], deg=1, cov=True)
    [mY_myFit, bY_myFit], covY_myFit = np.polyfit(y_dict[f'days_Set{my_set}'], y_dict[f'emittance_Set{my_set}_myFit'], deg=1, cov=True)
    # error on the slop of the fit
    errX_myFit = np.sqrt(np.diag(covX_myFit))[0]
    errX_myFit_Dx = np.sqrt(np.diag(covX_myFit_Dx))[0]
    errY_myFit = np.sqrt(np.diag(covY_myFit))[0]


    #print ("Slope X: " + str(mX))
    #print ("Intercept : " + str(bX))

    #print ("Slope Y: " + str(mY))
    #print ("Intercept : " + str(bY))

    # x-plane
    fig, ax = plt.subplots(1)

    xfmt = md.DateFormatter('%H:%M:%S')

    ax.plot(x_dict[f'days_Set{my_set}'], np.array(x_dict[f'emittance_Set{my_set}'])*1e6, 'o', c='b')
    ax.plot(x_dict[f'days_Set{my_set}'], np.array(x_dict[f'emittance_Set{my_set}_myFit'],)*1e6, 'o', c='C0')
    ax.plot(x_dict[f'days_Set{my_set}'], np.array(x_dict[f'emittance_Set{my_set}_myFit_with_Dx'],)*1e6, 'o', c='c')

    ax.plot(x_dict[f'days_Set{my_set}'], (np.array(x_dict[f'days_Set{my_set}'])*mX+bX)*1e6, c='b', label= 'parquet files: '+r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX*1e6/24:.2f}'+r'$\pm$'+f'{errX*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    ax.plot(x_dict[f'days_Set{my_set}'], (np.array(x_dict[f'days_Set{my_set}'])*mX_myFit+bX_myFit)*1e6, c='C0', label='myFit: '+r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX_myFit*1e6/24:.2f}'+r'$\pm$'+f'{errX_myFit*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    ax.plot(x_dict[f'days_Set{my_set}'], (np.array(x_dict[f'days_Set{my_set}'])*mX_myFit_Dx+bX_myFit_Dx)*1e6, c='c', label='myFit, Dx included: '+r'$\mathrm{d\epsilon_x/dt}$'+f'= {mX_myFit_Dx*1e6/24:.2f}'+r'$\pm$'+f'{errX_myFit_Dx*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    
    
    plt.ylim(1.00, 2.9)
    
    ax.set_title(f'Set {my_set}')
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$\mathrm{\epsilon_x \ [\mu m]}$')
    ax.legend(loc=4)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.gca().xaxis.set_major_formatter(xfmt)
    plt.grid(ls='--')
    plt.tight_layout()

    plt.savefig(f'./figures/myFit_11Oct2021/emitX_vs_time_Set{my_set}_compareFits.png', bbox_inches='tight')
    plt.close()

    # y-plane
    fig, ax = plt.subplots(1)

    xfmt = md.DateFormatter('%H:%M:%S')

    ax.plot(y_dict[f'days_Set{my_set}'], np.array(y_dict[f'emittance_Set{my_set}'])*1e6, 'o', c='r')
    ax.plot(y_dict[f'days_Set{my_set}'], np.array(y_dict[f'emittance_Set{my_set}_myFit'],)*1e6, '^', c='C1')
    
    ax.plot(y_dict[f'days_Set{my_set}'], (np.array(y_dict[f'days_Set{my_set}'])*mY+bY)*1e6, c='r', label='parquet files: '+r'$\mathrm{d\epsilon_y/dt}$'+f'= {mY*1e6/24:.2f}'+r'$\pm$'+f'{errY*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    ax.plot(y_dict[f'days_Set{my_set}'], (np.array(y_dict[f'days_Set{my_set}'])*mY_myFit+bY_myFit)*1e6, ls='--', c='C1', label='myFit: '+r'$\mathrm{d\epsilon_y/dt}$'+f'= {mY_myFit*1e6/24:.2f}'+r'$\pm$'+f'{errY_myFit*1e6/24:.2f} '+r'$\mathrm{[\mu m/h]}$')
    
    
    plt.ylim(1.00, 2.9)
    
    ax.set_title(f'Set {my_set}')
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$\mathrm{\epsilon_y \ [\mu m]}$')
    ax.legend(loc=4)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.gca().xaxis.set_major_formatter(xfmt)
    plt.grid(ls='--')
    plt.tight_layout()

    plt.savefig(f'./figures/myFit_11Oct2021/emitY_vs_time_Set{my_set}_compareFits.png', bbox_inches='tight')