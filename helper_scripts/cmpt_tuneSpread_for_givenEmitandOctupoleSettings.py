import numpy as numpy
from utils.my_functions import *
from utils.cmpt_TuneSpreads import *

Eb = 200e9  # eV
beta_0 = 0.999999  # cmpt it from the rest
gamma_0 = cmpt_gamma_0(Eb) # cmt it
betagamma = beta_0*gamma_0
print(f'betagamma = {betagamma}')

# Initial emittances
ey_n, ex_n = 3.339e-6, 5.366e-6 # m
ey_geom, ex_geom = ey_n/betagamma, ex_n/betagamma

# octupole settings
# detuning coefficients, [1/m]
# If you want to give it in k --> you need to connect it with the octupole matching.

a_xx = 0
a_xy  =  4509.24 #-1510.465536
a_yy =  -10000 #1080.331794 

Jy_rms, Jx_rms = ey_geom, ex_geom
#print(f'Initial rms actions: Jy ={Jy_rms}, Jx = {Jx_rms}')



Dq_rms = rms_amplitude_detuning_y_new(Jx_rms, Jy_rms, a_yy, a_xy)
print(f'The rms tune spread for ey={ey_n} m, ex={ex_n} m, axx={a_xx} 1/m, axy={a_xy} 1/m, axy={a_yy} 1/m is: \n {Dq_rms}')



