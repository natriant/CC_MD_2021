import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


if __name__ == '__main__':

    # get my current path
    my_path = os.getcwd().replace('\\', '/')

    # .h5 filename
    filename = '2021_10_13_CrabCavities.h5'
    filepath = os.path.join(my_path, filename).replace('\\', '/')

    # list of coast available
    coast_list = ['coast_01', 'coast_02', 'coast_03']
    coast_to_plot = 0

    # get data collected during the selected coast
    df = pd.read_hdf(filepath, coast_list[coast_to_plot])

    # Some plots
    plot_data = True
    save_plot = True

    if plot_data:

        fig, ax = plt.subplots(nrows=2, ncols=1, sharex='col')

        ax[0].errorbar(
            x=df['time'], y=df['tau_mean'], yerr=df['tau_std'],
            marker='.', markerfacecolor='C1', markeredgecolor='k',
            linestyle='', color='k',
            solid_capstyle='projecting', capsize=2.5,
            label=r'$\tau_{mean\ &\ std}$'
        )
        ax[0].plot(df['time'], df['tau_min'], marker='_', linestyle='', color='blue', alpha=0.5,
                   label=r'$\tau_{min}$')
        ax[0].plot(df['time'], df['tau_max'], marker='_', linestyle='', color='red', alpha=0.5,
                   label=r'$\tau_{max}$')

        myFmt = mdates.DateFormatter('%H:%M:%S')
        ax[0].xaxis.set_major_formatter(myFmt)
        # x-lim, beam debunched at ~ 15:45
        #ax[0].set_xlim(df['time'][1], df['time'][58])
        ax[0].set_ylim(1.8e-9, 2.3e-9)



        plt.xticks(rotation=30)

        ax[0].set_title('Bunch length')
        ax[0].set_ylabel('Bunch length [s]')
        fig.suptitle('2021.10.13 - Crab Cavities tests: ' + str(coast_list[coast_to_plot]), fontsize=16)

        ax[0].grid(linestyle="--")
        ax[0].legend(loc=2)

        ax[1].errorbar(
            x=df['time'], y=df['peak_mean'], yerr=df['peak_std'],
            marker='.', markerfacecolor='C1', markeredgecolor='k',
            linestyle='', color='k',
            solid_capstyle='projecting', capsize=2.5,
            label=r'$Peak_{mean\ &\ std}$'
        )
        ax[1].plot(df['time'], df['peak_min'], marker='_', linestyle='', color='blue', alpha=0.5,
                   label=r'$Peak_{min}$')
        ax[1].plot(df['time'], df['peak_max'], marker='_', linestyle='', color='red', alpha=0.5,
                   label=r'$Peak_{max}$')

        myFmt = mdates.DateFormatter('%H:%M:%S')
        ax[1].xaxis.set_major_formatter(myFmt)
        # x-lim, beam debunched at ~ 15:45
        #ax[1].set_xlim(df['time'][1], df['time'][len(df['time'])])
        ax[1].set_ylim(0.2, 0.4)

        plt.xticks(rotation=30)

        ax[1].set_title('Peak amplitude')
        ax[1].set_xlabel('Acquisition time')
        ax[1].set_ylabel('Voltage [a.u.]')

        ax[1].grid(linestyle="--")
        ax[1].legend(loc=3)

        plt.tight_layout()

        if save_plot:
            fig_name = str(coast_list[coast_to_plot]) + '.png'
            plt.savefig(fig_name)

        plt.show()
