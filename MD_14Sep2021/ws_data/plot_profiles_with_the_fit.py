'''
1. Open the parquet files.
2. Compute the emittances with a gaussian fit on the profiles.
3. Plot the emittance evolution along with a polynomial fit of 1st degree.
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
          'figure.figsize': (10, 7),
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
my_variables = ['emittance_Set', 'profiles_Set', 'positions_Set', 'betagamma_Set', 'timestamp_Set']
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
                x_dict[f'timestamp_Set{my_set}'].append(data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6)  # Convert UNIX time to days since Matplotlib epoch.
                for variable in my_variables[:-2]:
                    x_dict[f'{variable}{my_set}'].append(data[f'{variable}{my_set}'][0][0])
                   
            if '41677.V' in filename:
                y_dict[f'betagamma_Set{my_set}'].append(data[f'betagamma_Set{my_set}'])
                y_dict[f'timestamp_Set{my_set}'].append(data['cycleStamp'][0]/1e9+ data[f'acq_time_Set{my_set}'][0]/1e6)
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

    y_dict[f'Set{my_set}_myFitresults'] = []
    x_dict[f'Set{my_set}_myFitresults'] = []
    
    for i in range(len(y_dict[f'positions_Set{my_set}'])): # iterate over the acquisitions for each setting, y-plane
        pop_t_y, errors = fitGauss(y_dict[f'positions_Set{my_set}'][i], y_dict[f'profiles_Set{my_set}'][i], fit_func='Gauss5p')
        A, mu, sigma, offset = pop_t_y[0], pop_t_y[1], pop_t_y[2], pop_t_y[3]
        ey = getEmittance(sigma, beta['SPS.BWS.41677.V'], y_dict[f'betagamma_Set{my_set}'][i]) # returns the emittance in [um]
        err= errors[-1]/sigma # from MD 2018. cernbox/CC_MDs_2018/27July2020/MD5_5Sep2018_beam_profiles_vs_emitBU_coast2.ipynb. Added by Natalia 14.10.2021
        
        y_dict[f'emittance_Set{my_set}_myFit'].append([ey[0], err])
        results = [A, mu, sigma, offset]
        y_dict[f'Set{my_set}_myFitresults'].append(results)

    for i in range(len(x_dict[f'positions_Set{my_set}'])): # iterate over the acquisitions for each setting, x-plane
        pop_t_x, errors = fitGauss(x_dict[f'positions_Set{my_set}'][i], x_dict[f'profiles_Set{my_set}'][i], fit_func='Gauss5p')
        A, mu, sigma, offset = pop_t_x[0], pop_t_x[1], pop_t_x[2], pop_t_x[3]
        ex = getEmittance(sigma, beta['SPS.BWS.51637.H'], x_dict[f'betagamma_Set{my_set}'][i])
        ex_Dx = getEmittance_with_Dx(pop_t_x[2],  beta['SPS.BWS.51637.H'], x_dict['betagamma_Set1'][0], Dx['SPS.BWS.51637.H'], dpp) # returns the emittance in [um]
        err= errors[-1]/sigma
        
        x_dict[f'emittance_Set{my_set}_myFit'].append([ex[0], err])
        x_dict[f'emittance_Set{my_set}_myFit_with_Dx'].append([ex_Dx[0], err])
        results = [A, mu, sigma, offset]
        x_dict[f'Set{my_set}_myFitresults'].append(results)


def yRefFun(x, my_A, my_mu, my_sigma, my_offset):
    return my_A * np.exp(-0.5 * (x - my_mu) ** 2 / my_sigma ** 2) + my_offset

for i in range(len(x_dict[f'positions_Set1'])): # iterate over the acquisitions for each setting
    fig, ax = plt.subplots()
    
    for my_set in range(1,3):

        results = y_dict[f'Set{my_set}_myFitresults'][i]
        A, mu, sigma, offset = results[0], results[1], results[2], results[3]

        xRef = y_dict[f'positions_Set{my_set}'][i]-mu
        yRef = y_dict[f'profiles_Set{my_set}'][i]-offset

        k = yRefFun(xRef, A, 0, sigma, 0)
        max_to_normalize = max(k)
       
        ey, err_y = y_dict[f'emittance_Set{my_set}_myFit'][i][0], y_dict[f'emittance_Set{my_set}_myFit'][i][1]
        ax.plot(xRef, yRef, 'o', ms=4, C=f'C{my_set-1}', label=r'$\mathrm{\epsilon_y}$'+f',Set {my_set}={ey:.3f}' +r'$\pm /$'+ f'{err_y:.3f}'+ r'$\mathrm{[\mu m}]$')
        ax.plot(xRef, k, C=f'C{my_set-1}')
        ax.set_xlim(-5, 5)
        ax.set_ylim(-10, 4000)

        ax.legend(loc=2)
    title = str(datetime.datetime.fromtimestamp(x_dict[f'timestamp_Set{my_set}'][i]))[:19]
    ax.set_title(title)
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'41677.V_IN_OUT_{title[10:19]}.png', bbox_inches='tight')