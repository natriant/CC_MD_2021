'''
Measured growth: sum x + y planes. Average IN and OUT scan
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

params = {'legend.fontsize': 20,
          'figure.figsize': (9.5, 8.5),
          'axes.labelsize': 27,
          'axes.titlesize': 23,
          'xtick.labelsize': 27,
          'ytick.labelsize': 27,
          'image.cmap': 'jet',
          'lines.linewidth': 1,
          'lines.markersize': 10,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)

measured_2018 =  [18.07518005914758, 0.4716569640922441, 1.1793448259787538, 1.3245910844084205, 2.364473683403215, 9.37461985936439]
measured_2018_err =[1.0258202495920914, 0.17445331926616509, 0.14799288858098972, 0.2728254877629634, 0.1891975770128703, 0.6489717538012192]
expected_2018 = [40.04166677228909, 1.8364995968548514, 5.291054200334389, 5.332086619494902, 14.883216955922814, 48.0104746749581]

back_x, back_y = 1.09, 0.72

measured_2021 = np.array([3.425, 3.57, 7.285, 4.455])-(back_x+back_y)
measured_2021_err = [0.81, 0.65, 1.17, 0.82]

df_expected_2021 = pd.read_pickle('../output/expected_growths_CC_MD_13Oct2021_Vcc1_1MV.pkl')

coast1_setting1 = df_expected_2021['expected dey/dt AN [um/h]'][0]+ df_expected_2021['expected dey/dt PN [um/h]'][0]
coast1_setting2 = df_expected_2021['expected dey/dt AN [um/h]'][1]+ df_expected_2021['expected dey/dt PN [um/h]'][1]
coast1_setting3 = df_expected_2021['expected dey/dt AN [um/h]'][2]+ df_expected_2021['expected dey/dt PN [um/h]'][2]
coast1_setting3b = df_expected_2021['expected dey/dt AN [um/h]'][3]+ df_expected_2021['expected dey/dt PN [um/h]'][3]
coast2_setting1 = df_expected_2021['expected dey/dt AN [um/h]'][4]+ df_expected_2021['expected dey/dt PN [um/h]'][4]
coast3_setting1 = df_expected_2021['expected dey/dt AN [um/h]'][5]+ df_expected_2021['expected dey/dt PN [um/h]'][5]
coast3_setting2 = df_expected_2021['expected dey/dt AN [um/h]'][6]+ df_expected_2021['expected dey/dt PN [um/h]'][6]
coast3_setting3 = df_expected_2021['expected dey/dt AN [um/h]'][7]+ df_expected_2021['expected dey/dt PN [um/h]'][7]



#expected_2021 = [2, 5, 10, 5]
expected_2021=[coast1_setting1, coast1_setting2, coast1_setting3, coast2_setting1]

#expected_2021_oct = [10, 5, 10, 10]
expected_2021_oct = [coast1_setting3b, coast3_setting1, coast3_setting2, coast3_setting3]
measured_2021_oct = np.array([11.06, 4.165, 7.93, 7.41])-(back_x+back_y)
measured_2021_oct_err =[2.47, 0.64, 0.99, 1.43]


fig, ax = plt.subplots(figsize=(10, 8))
ax.errorbar(expected_2018, measured_2018, yerr= measured_2018_err, marker='o', ls='', markersize=13, capsize=5, label='2018')
ax.errorbar(expected_2021, measured_2021, yerr= measured_2021_err, marker='o', ls='', markersize=13, capsize=5, label='2021')
ax.errorbar(expected_2021_oct, measured_2021_oct, yerr= measured_2021_oct_err, marker='o', ls='', markersize=13, capsize=5, label='2021, LOD')

ax.grid(ls='--')
ax.set_xlabel('Theoretically expected growth [um/h]')
ax.set_ylabel('Measured growth [um/h]')
ax.set_xlim(-2, 50)
ax.set_ylim(-2, 50)


ax.set_title('Preliminary', fontsize=27)

ax.plot([-50, 50], [-50, 50], ls="--", c=".3", lw=3)
ax.plot([-100, 100], [-50, 50], ls="--", c="moccasin", lw=3)#, label=r'$\mathrm{y=\frac{1}{2}x}$')
ax.plot([-100, 100], [-25, 25], ls="--", c="skyblue", lw=3)#, label=r'$\mathrm{y=\frac{1}{4}x}$')

ax.set_yticks(np.arange(0,60,10))
ax.set_xticks(np.arange(0,60,10))

#ax.set_xlim(-2, 20)
#ax.set_ylim(-2, 20)


ax.set_aspect('equal', adjustable='box')
plt.legend(loc=2, fontsize=27)
plt.tight_layout()
#plt.axis('equal')
plt.savefig('CC_MD_emitGrowth_2018vs2021_LOD_v3.png',  bbox_inches="tight")
