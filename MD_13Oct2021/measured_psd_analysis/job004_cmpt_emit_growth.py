import pickle as pkl
import pandas as pd
import sys

sys.path.append('../../utils/')
from cmptTheoreticalEmitGrowth import *
from bunchLengthConversions import *
from coordinatesConversions import *

clight = 299792458 # m/s

# Machine parameters
betay = 73.81671646 # 73 m at CC2, 76 m at CC1
Vcc = 1e6 # V
f_CC = 400e6 # CC frequency in Hz
circumference = 6911.5623 # m

# Beam parameters
Eb = 200e9 # 270e9 eV
gamma_0 = 213.16 # for 200 GeV  # 287.8 for 270 GeV
beta_0 = np.sqrt(1 - 1/gamma_0**2)
frev = 299792458/circumference # Hz
tau = 2.0e-9 # 4 sigma_t [s] # 1.7e-9
sigma_z = clight*tau/4 #0.155  # m
print(f'sigma_z = {sigma_z} m')


df =  pd.read_pickle('psd_at_fb_CC_MD_13Oct2021.pkl')
df_keys = list(df.keys())

for i in range(len(df)): # enumerate over the rows of the data frame
	name = df.iloc[i].name # coast setting 
	AN_psd = df.iloc[i][0][0] #[0][0] 0--> AN, 0 --> mean, no std, in 1/Hz
	PN_psd = df.iloc[i][1][0] #[0][0] 1--> PN, 0 --> mean, no std, in rad^2/Hz

	# Compute the correction factor due to the bunch length
	sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_CC)
	CDeltaphi = cmpt_bunch_length_correction_factor(sigma_phi, 'PN')
	CDeltaA = cmpt_bunch_length_correction_factor(sigma_phi, 'AN')


	expected_growth_PN = emit_growth_phase_noise(betay, Vcc, frev, Eb, CDeltaphi, PN_psd, one_sided_psd=True)*beta_0*gamma_0 # m/s
	expected_growth_AN = emit_growth_amplitude_noise(betay, Vcc, frev, Eb, CDeltaA, AN_psd, one_sided_psd=True)*beta_0*gamma_0 # m/s
	
	print(name)
	print(f'PN dey/dt {expected_growth_PN*3600*1e6:.3f} [um/h]')
	print(f'AN dey/dt {expected_growth_AN*3600*1e6:.3f} [um/h]')