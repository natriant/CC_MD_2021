import numpy as np

def chromatic_tune_spread(dpp_rms, order, Qp):
    if order == 1:
        print('1st chromaticity')
        dqy = Qp * dpp_rms
    if order == 2:
        print('2nd order chromaticity')
        dqy = Qp * (dpp_rms ** 2)
    return dqy


def amplitude_detuning_x(Jx, Jy, a_xx, a_xy):
        # compatible with pyheadtail
        return a_xx * 2 * Jx + a_xy * 2 * Jy


def amplitude_detuning_y(Jx, Jy, a_yy, a_xy):
    return a_yy * 2 * Jy + a_xy * 2 * Jx


def rms_amplitude_detuning_x(Jx, Jy, a_xx, a_xy):
        # compatible with pyheadtail
        # Jx, Jy: arrays
        return 2*np.sqrt((a_xx*np.std(Jx))**2+(a_xy*np.std(Jy))**2)


def rms_amplitude_detuning_y(Jx, Jy, a_yy, a_xy):
    # Jx, Jy: arrays
    return 2*np.sqrt((a_yy*np.std(Jy))**2+(a_xy*np.std(Jx))**2)



def rms_amplitude_detuning_x_new(Jx_rms, Jy_rms, a_xx, a_xy):
        # compatible with pyheadtail
        # Jx, Jy: arrays
        return 2*np.sqrt((a_xx*Jx_rms)**2+(a_xy*Jy_rms)**2)

def rms_amplitude_detuning_y_new(Jx_rms, Jy_rms, a_yy, a_xy):
    # Jx, Jy: arrays
    return 2*np.sqrt((a_yy*Jy_rms)**2+(a_xy*Jx_rms)**2)

