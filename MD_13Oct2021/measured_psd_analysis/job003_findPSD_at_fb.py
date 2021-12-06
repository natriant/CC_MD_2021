import numpy as np
import pandas as pd
import pickle as pkl
import sys
from utils import *
import matplotlib.pyplot as plt

saveflag=False # True to save data

#### SPS parameters ######
circumference = 6911.5623
frev = 299792458/circumference
#################################3

# Load the measured spectrum in a data frame
path_to_data = './measured_noise_files/2021/'

coasts = [1,2,3] # given number of coast. choose 1 2 or 3

settings = [1,2,3] #[1, 2, 3] # number of settings for this coast
noise_types = ['AM', 'PM']

column_names = ['AM [1/Hz]', 'PM [rad^2/Hz]']
df_results = pd.DataFrame(columns = column_names, dtype=object)
index_list = []

for coast in coasts:
    for setting in settings:
        to_append = []
        p1=0
        for noise_type in noise_types:
            file_name = f'COAST{coast}-SETTING{setting}-CC1-{noise_type}'	
            print(file_name)
            try:
                df = pd.read_csv(path_to_data + file_name+'.csv')
                dict_keys = df.keys()       
            

                L_log_load = np.array(df[dict_keys[1]]) # PSD in dBc/Hz
                freq_load = np.array(df[dict_keys[0]])  # equally spaced in logarithmic scale

                # A. Convert SSB, 10Log_10L [dBc/Hz] to G_yy [rad^2/Hz]
                Gyy = ssb_2_dsb(np.array(L_log_load))  # PSD in rad^2/Hz

                N = int(1e5)+1 #100001 #500001  # number of samples --> number of simulated turns

                t = np.linspace(0, N / frev, N)  # time in seconds
                delta_t = t[1] - t[0]  # this should be fixed
                freq = np.fft.fftfreq(N, delta_t)  # positive and negative frequencies
                delta_f = freq[1] - freq[0]
                print('Delta t = {} Hz'.format(delta_t))
                print('Delta f = {} Hz'.format(delta_f))

                # for the linear interpolation, we need only the positive frequencies
                freq_pos = freq[1:N // 2+1] # https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html
                Gyy_linear_interpolation = np.interp(freq_pos, freq_load, Gyy)

                ### Create the positive spectral components of the two-sides spectrum S_yy
                Syy = Gyy_linear_interpolation/2
                
                fb, step = 0.18*frev, 500
                my_mean, my_std = find_average_psd_around_freq(freq, Syy, fb, step)
                print(f'{noise_type} noise spectrum: The mean PSD value around {fb} Hz over a range of {step} Hz is {my_mean} +-{my_std} rad^2/Hz')
                to_append.append([my_mean, my_std])


                # fill in index list
                if noise_type=='AM':
                    index_list.append(f'coast{coast}-setting{setting}')
                
            except OSError as err:
                p1=1
                print("OS error: {0}".format(err))
                pass
        if p1==0:
            df_length = len(df_results)
            df_results.loc[df_length] = to_append

df_results.index = index_list

print(df_results)
if saveflag:
    df_results.to_pickle('./psd_at_fb_CC_MD_13Oct2021.pkl')