import numpy as np


def compute_PSD_whiteNoise(N, frev, muNoise, stdNoise):
    t = np.linspace(0, N / frev, N)  # time in seconds
    delta_t = t[1] - t[0]
    freq = np.fft.fftfreq(N, delta_t)
    delta_f = freq[1] - freq[0]
    print('sampling Delta t = {} s'.format(delta_t))
    print('frequency resolution, Delta f = {} s'.format(delta_f))

    #### Iterate over 1000 times to increase the precision of the FFT
    fft_list = []
    for i in range(1000):
        y_noise = np.random.normal(muNoise, stdNoise, N)
        fft = np.fft.fft(y_noise)
        fft_list.append(fft)
    mean_dft = np.mean(np.abs(fft_list)**2, axis=0)
    psd = mean_dft/(delta_f*(N**2))

    return psd, freq


def compute_PSD(signal, N, frev):
    t = np.linspace(0, N / frev, N)  # time in seconds
    delta_t = t[1] - t[0]
    freq = np.fft.fftfreq(N, delta_t)
    delta_f = freq[1] - freq[0]
    print('sampling Delta t = {} s'.format(delta_t))
    print('resolution frequency, Delta f = {} s'.format(delta_f))
    fft = np.fft.fft(signal)
    psd = np.abs(fft)**2/(delta_f*(N**2))

    return psd, freq


def ssb_2_dsb(L):
    """
    Convert the single sideband, SBB, measurement of the noise, L(f) in dBc/Hz, to the double
    sideband density, DSB,  S(f) in rad^2/Hz. According to the IEEE standards:
    S(f) = 2*10^(L(f)/10) in rad^2/Hz, where L(f) in dBc/Hz.
    IMPORTANT: In the definitions here, both L(f) and S(f) are one sided, ie only positive frequencies.
    """

    S = 2*10**(L/10)
    return S


def find_closest_to_given_value(my_array, my_value):
    k = my_array[min(range(len(my_array)), key=lambda i: abs(my_array[i] - my_value))]
    return k


def find_index_for_value(my_array, my_value):
    my_index = [i for i in range(len(my_array)) if my_array[i] == my_value]
    return my_index[0]


def find_average_psd_around_freq(freq_list, sxx_list, my_freq, step):
    closest_to_fb = find_closest_to_given_value(freq_list, my_freq)
    closest_to_fb_index = find_index_for_value(freq_list, closest_to_fb)

    # Define the range of frequency values for averaging the PSD values
    fb_min, fb_max = my_freq - step, my_freq + step
    closest_to_fb_min = find_closest_to_given_value(freq_list, fb_min)
    closest_to_fb_min_index = find_index_for_value(freq_list, closest_to_fb_min)
    closest_to_fb_max = find_closest_to_given_value(freq_list, fb_max)
    closest_to_fb_max_index = find_index_for_value(freq_list, closest_to_fb_max)

    my_mean = np.mean(sxx_list[closest_to_fb_min_index:closest_to_fb_max_index])
    my_std = np.std(sxx_list[closest_to_fb_min_index:closest_to_fb_max_index])

    return my_mean, my_std
