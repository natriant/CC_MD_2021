import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (10, 7),
          'axes.labelsize': 25,
          'axes.titlesize': 23,
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


if __name__ == '__main__':

    # get my current path
    my_path = os.getcwd().replace('\\', '/')

    # .h5 filename
    filename = '2021_10_13_CrabCavities.h5'
    filepath = os.path.join(my_path, filename).replace('\\', '/')

    # list of coast available
    coast_list = ['coast_01', 'coast_02', 'coast_03']
    coast_to_plot = 2

    # get data collected during the selected coast
    df = pd.read_hdf(filepath, coast_list[coast_to_plot])

    # Some plots
    plot_data = True
    save_plot = True

    if plot_data:

        #fig, ax = plt.subplots(nrows=2, ncols=1, sharex='col')
        fig, ax = plt.subplots(nrows=1, ncols=1)#, sharex='col')

        ax.errorbar(
            x=df['time'], y=df['tau_mean']*1e9, yerr=df['tau_std']*1e9,
            marker='.', markerfacecolor='C1', markeredgecolor='k',
            linestyle='', color='k',
            solid_capstyle='projecting', capsize=2.5,
            label=r'$\tau_{mean\ &\ std}$'
        )
        ax.plot(df['time'], df['tau_min']*1e9, marker='_', linestyle='', color='blue', alpha=0.5,
                   label=r'$\tau_{min}$')
        ax.plot(df['time'], df['tau_max']*1e9, marker='_', linestyle='', color='red', alpha=0.5,
                   label=r'$\tau_{max}$')

        myFmt = mdates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(myFmt)
        # x-lim, beam debunched at ~ 15:45
        #ax[0].set_xlim(df['time'][1], df['time'][58])
        ax.set_ylim(1.8, 2.3)



        plt.xticks(rotation=30)

        ax.set_title('2021.10.13 - Crab Cavities tests: ' + str(coast_list[coast_to_plot]))
        ax.set_ylabel('4σ'+r'$\mathrm{_t}$'+' [ns]')
        ax.set_xlabel('Time')
        #fig.suptitle('2021.10.13 - Crab Cavities tests: ' + str(coast_list[coast_to_plot]), fontsize=21)

        ax.grid(linestyle="--")
        ax.legend(loc=2)
        
        
        plt.tight_layout()

        if save_plot:
            fig_name = str(coast_list[coast_to_plot]) + 'v2.png'
            plt.savefig(fig_name)

        plt.show()
