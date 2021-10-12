# Compute the expected emittance growth from the theoretical model of Themis and Philippe.


import numpy as np
from utils.my_functions import *

### Beam and machine parameters
betay = 73.82  # m in CC2, ~76 in CC1 (MAD-X)
Vcc = 1e6  # V
frev = 43.45e3  # Hz
Eb = 270e9  # eV
beta_0 = 0.999999  # cmpt it from the rest
gamma_0 = cmpt_gamma_0(Eb) # cmt it
clight = 299792458  # light speed in meters/second
f_RF = 400.789e6  # CC frequency in Hz
f_CC = 200e6

tau = 2.2e-9 # 4 sigma_t [s] # 1.7e-9
sigma_z = clight*tau/4 #0.155  # m
print(f'sigma_z = {sigma_z} m')

PN = -114.8 # dbc/Hz
AN = 0

##############################


PSD_phi, PSD_A = ssb_2_dsb(PN), 0#ssb_2_dsb(AN) 
print(f'Phase noise in rad^2/Hz: {PSD_phi}')
print(f'Amplitude noise in rad^2/Hz: {PSD_A}')


# Compute the correction factor due to the bunch length
sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_CC)
CDeltaphi = cmpt_bunch_length_correction_factor(sigma_phi, 'PN')
CDeltaA = cmpt_bunch_length_correction_factor(sigma_phi, 'AN')
print(f'CDeltaPhi = {CDeltaphi}, CDeltaA = {CDeltaA}')

expected_growth_PN = emit_growth_phase_noise(betay, Vcc, frev, Eb, CDeltaphi, PSD_phi, one_sided_psd=True)*beta_0*gamma_0 # m/s
expected_growth_AN = emit_growth_amplitude_noise(betay, Vcc, frev, Eb, CDeltaA, PSD_A, one_sided_psd=True)*beta_0*gamma_0 # m/s
expected_growth = expected_growth_PN + expected_growth_AN


conversion_factor = 1e-3*3600 # convert to um/h
print(f'Expected growth from phase noise: {expected_growth_PN*1e9*conversion_factor} um/h')
print(f'Expected growth from amplitude noise: {expected_growth_AN*1e9*conversion_factor} um/h')
print(f'Expected growth total: {expected_growth*1e9*conversion_factor} um/h')




