import pickle as pkl
import pandas as pd
import sys

sys.path.append('../../../utils/')
from cmptTheoreticalEmitGrowth import *
from bunchLengthConversions import *
from coordinatesConversions import *

### Import measured bunch length and psd #####
taus_df = pd.read_pickle('mean_tau_summary_CC_MD_13Oct2021.pkl')
psd_df =  pd.read_pickle('psd_at_fb_CC_MD_13Oct2021.pkl')

# re-form psd data frame such as it is compatible with the taus data frame
line_1 = psd_df.iloc[2]
line_2 = psd_df.iloc[-1]
psd_df.loc[len(psd_df)] = line_1
psd_df.loc[len(psd_df)] = line_2

index_list = ['coast1-setting1', 'coast1-setting2', 'coast1-setting3',  'coast1-setting3b', 'coast2-setting1', 'coast3-setting1', 'coast3-setting2', 'coast3-setting3']

psd_df.index = index_list

a, b, c, d, e = psd_df.iloc[3].copy(), psd_df.iloc[4].copy(),  psd_df.iloc[5].copy(),  psd_df.iloc[6].copy(),  psd_df.iloc[7].copy() 
psd_df.iloc[3], psd_df.iloc[4], psd_df.iloc[5], psd_df.iloc[6], psd_df.iloc[7] = d, a, b, c, e


print(taus_df)
print(psd_df)



### Machine and other beam parameters

clight = 299792458 # m/s

# Machine parameters
betay = 76.07 #73.81671646 # 73 m at CC2, 76 m at CC1
Vcc = 1e6 # V
f_CC = 400e6 # CC frequency in Hz
circumference = 6911.5623 # m

# Beam parameters
Eb = 200e9 # 270e9 eV
gamma_0 = 213.16 # for 200 GeV  # 287.8 for 270 GeV
beta_0 = np.sqrt(1 - 1/gamma_0**2)
frev = 299792458/circumference # Hz


# Save theoretically expected growth rates
column_names = ['expected dey/dt AN [um/h]', 'expected dey/dt PN [um/h]'] 
df_results = pd.DataFrame(columns = column_names, dtype=object)

for i in range(len(psd_df)):
	name = psd_df.iloc[i].name # coast setting
	
	AN_psd = psd_df['AM [1/Hz]'][i][0] # --> mean, no std, in 1/Hz
	PN_psd = psd_df['PM [rad^2/Hz]'][i][0] # --> mean, no std, in rad^2/Hz


	tau = taus_df['tau [ns], average'][i+1]
	#tau = taus_df['tau [ns], initial'][i+1] # 
	sigma_z = clight*tau/4 #0.155  # m
	print(f'sigma_z = {sigma_z} m')


	# Compute the correction factor due to the bunch length
	sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_CC)
	CDeltaphi = cmpt_bunch_length_correction_factor(sigma_phi, 'PN')
	CDeltaA = cmpt_bunch_length_correction_factor(sigma_phi, 'AN')


	expected_growth_PN = emit_growth_phase_noise(betay, Vcc, frev, Eb, CDeltaphi, PN_psd, one_sided_psd=True)*beta_0*gamma_0 # m/s
	expected_growth_AN = emit_growth_amplitude_noise(betay, Vcc, frev, Eb, CDeltaA, AN_psd, one_sided_psd=True)*beta_0*gamma_0 # m/s
	
	print(name)
	print(f'PN dey/dt {expected_growth_PN*3600*1e6:.3f} [um/h]')
	print(f'AN dey/dt {expected_growth_AN*3600*1e6:.3f} [um/h]')

	df_results.loc[i] = [expected_growth_AN*3600*1e6, expected_growth_PN*3600*1e6]

df_results.index = index_list
print(df_results)

saveflag=False
if saveflag:
    df_results.to_pickle('./output/expected_growths_CC_MD_13Oct2021.pkl')

