from utils.MR import *
from datetime import datetime
import glob
import sdds
import matplotlib.pyplot as plt

# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (8, 7),
          'axes.labelsize': 25,
          'axes.titlesize': 21,
          'xtick.labelsize': 23,
          'ytick.labelsize': 23,
          'image.cmap': 'jet',
          'lines.linewidth': 3,
          'lines.markersize': 10,
          'font.family': 'sans-serif'}


plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)
plt.ion()

def get_average_profile(my_profiles_amp): #  from utils/MR.py
    average_profile = np.mean(my_profiles_amp, axis=0) 
    baseline = np.mean(average_profile[0:20]) 
    average_profile -= baseline
    return average_profile

path2files='SMR.SCOPE13.CH01@Acquisition/' # path to MR data, they should be unziped (.sdds)
files_list = sorted(glob.glob(path2files+'*')) # create a list with all the MR data

mr_overview = pickle.load(open('mr_overview.pkl', 'rb'))
sigma_t = mr_overview['LHCMD4']['sigma_t'][120:287]
mu =  mr_overview['LHCMD4']['mu0'][120:287]
amp = mr_overview['LHCMD4']['ampl'][120:287]
# [120:290] Acquisitions that correspond to the time of the emittance growth measurements. [120:287] chosen to analyse exactly 1h of measurements.
# Index : 287 to be exactly 1h of measurements (real time up to 290)

for i, my_file in enumerate(files_list[120:287]):

    sdds_data = sdds.read(my_file)
    profiles_amp =  sdds_data.values['value']
    triggerStamps = sdds_data.values['triggerStamp']
    start_str = str(datetime.fromtimestamp(triggerStamps [10]/1e9))
    end_str = str(datetime.fromtimestamp(triggerStamps [-1]/1e9))
    
    average_profile = get_average_profile(profiles_amp)
    
    time_axis = np.float_(range(profiles_amp.shape[1]))*sdds_data.values['sampleInterval']-profiles_amp.shape[1]*sdds_data.values['sampleInterval']/2.

    # Get fit
    p=[amp[i], mu[i], sigma_t[i]]
    fitfunc = lambda p, x: p[0]*np.exp(-((x-p[1])**2)/(2.*p[2]**2)) # Gaussian function
    my_fit =fitfunc(p, time_axis)
    
    # Get difference between two datetimes in milliseconds
    date_format_str = '%Y-%m-%d %H:%M:%S.%f'
    start = datetime.strptime(start_str, date_format_str)
    end =   datetime.strptime(end_str, date_format_str)
    # Get the interval between two datetimes as timedelta object
    diff = end - start
    # Get the interval in milliseconds
    diff_in_milli_secs = diff.total_seconds() * 1000
    
    my_title = start_str[:19] # Y-M-D-H-M-S
    fig, ax = plt.subplots()
    ax.plot(time_axis, average_profile, c='g', label= f'average profile over {diff_in_milli_secs:.2f} [ms]')  
    ax.plot(time_axis, my_fit, c='k', lw=2, ls='--', label=f'gaussian fit,\n' + r'$\mathrm{4\sigma_t}$'+f'={4*sigma_t[i]:.2} [ns]')
   
    
    ax.set_xlim(-10, 10)
    ax.set_ylim(-20, 6000)
    
    ax.set_title(f'{my_title}')
    ax.set_ylabel('Amplitude [a.u.]') # Maybe in Volts
    ax.set_xlabel('Time [ns]') # Not entirely sure that is time in [ns]
    ax.legend(loc=2, frameon=False)
    ax.grid(ls='--')
    
    
    
    plt.tight_layout()
    plt.savefig(f'./figures/profiles/profiles_GaussFit_{my_title}.png', bbox_inches='tight')
    plt.close()