'''
From: link to swan
Link to overleaf: https://www.overleaf.com/2439138969dpmgdqdfhzny
Original script: Project_thesis/scripts_for_simple_figures_plots/measured_psd_to_time_series
'''

import os
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from utils import *

# List all the noise files
path_to_data = './measured_noise_files/2021/'
files_list = os.listdir(path_to_data)
print(files_list)
save_flag = True

files2ignore_list = ['COAST1-SETTING0-CC1-PM-NOEXCIT.csv', 'COAST1-SETTING0-CC1-AM.csv'] # ignore the files for no excitation
# to reduce the statistical uncertainty of the simulations, which is introduced from the random angles theta, 5 different sets from the
# same measured psd are created.
n_sets = 5 # int(sys.argv[1])  

# A. Load data
for n_set in range(0, n_sets):
    for file in files_list:
        if file not in files2ignore_list:
            df = pd.read_csv(path_to_data + file)
            dict_keys = df.keys()
            
            L_log_load = np.array(df[dict_keys[1]]) # PSD in dBc/Hz
            freq_load = np.array(df[dict_keys[0]])  # equally spaced in logarithmic scale

            # B. Generate time series
            # B1. Convert SSB, 10Log_10L [dBc/Hz] to G_yy [rad^2/Hz]
            Gyy = ssb_2_dsb(np.array(L_log_load))  # PSD in rad^2/Hz

            # B2. Linear interpolation of the frequency array, such as the values are equally spaced linearly

            clight = 299792458  # m/s
            circumference = 6911.5623
            frev = clight / circumference  # Hz
            print(frev)

            N = int(1e5)+1 #100001 #500001  # number of samples --> number of simulated turns

            t = np.linspace(0, N / frev, N)  # time in seconds
            delta_t = t[1] - t[0]  # this should be fixed
            # freq = np.linspace(0, N/t[-1], N) # [0, 2frev]
            freq = np.fft.fftfreq(N, delta_t)  # positive and negative frequencies
            print(len(freq))

            delta_f = freq[1] - freq[0]
            print('Delta t = {} Hz'.format(delta_t))
            print('Delta f = {} Hz'.format(delta_f))

            # for the linear interpolation, we need only the positive frequencies
            freq_pos = freq[1:N // 2+1] # https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html
            print(freq_pos[-1]) # tests if the last value is positive
            print('positive frequencies ={}'.format(len(freq_pos)))
            Gyy_linear_interpolation = np.interp(freq_pos, freq_load, Gyy)

            # B3. Create the positive spectral components of the two-sides spectrum S_yy
            Syy = Gyy_linear_interpolation/2

            # B4. Compute the amplitude of the spectral components of the Fourier transform, only for the postive part of the spectrum
            y_fft_amplitude = np.sqrt(Syy * delta_f * (N ** 2))

            # B5. Create random angles - for sets > 1 
            # The following lines ensure that you draw a random number using a different seed for each turn
            seed = 12345 + n_set
            rng = np.random.default_rng(seed)

            _ = rng.uniform(size=int(len(freq_pos)))
            theta = 2 * np.pi * rng.uniform(0, 1, int(len(freq_pos)))

            # B5. Create random angles - for 1 set
            #theta = 2 * np.pi * np.random.uniform(0, 1, int(len(freq_pos)))


            # B6. Create the new FFT
            y_fft_amplitude_angles = y_fft_amplitude * np.exp(1j * theta)  # note, this does contain the zero frequency term. a[1:N//2]

            # B7. Hermitian, to obtain real signal as a result
            fft_temp = np.concatenate((np.empty_like(y_fft_amplitude_angles), np.empty_like(y_fft_amplitude_angles)), axis=0) # constract an array with the spectral components in the order of the numpy.fft function. We need the concatenate and empty like function such as we create an array who can taki imaginary numbers as input
            fft_temp = np.insert(fft_temp, 0, 0, axis=0) # insert DC in the beggining of the array. We set the zero frequency term to be only real and equals zero.
            fft_temp[1:(N // 2 + 1)] = y_fft_amplitude_angles  # from 1 to 5000
            fft_temp[N // 2 + 1:] = np.conj(np.flip(y_fft_amplitude_angles))  # from 5001 to 10001
            len(fft_temp)

            signal_new = np.fft.ifft(fft_temp)  # time series, noise kicks
        
            if save_flag:
                with open(f'./output/time_series_from_{file}_set{n_set}.pkl', 'wb') as f2:
                    pkl.dump(np.real(signal_new), f2)
                f2.close()
        else:
            print(f'{file} ignored')
            