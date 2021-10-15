import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


# From: https://gitlab.cern.ch/scripting-tools/injectors-beam-monitoring/-/blob/master/injectors_beam_monitoring/analysis/wire_scanner_analysis.py
def fitGauss(pos, prof, pos_min=-np.inf, pos_max=np.inf, fit_func='Gauss4p'):
    assert fit_func in ["Gauss4p", "Gauss5p"]
    i = np.where((pos > pos_min) & (pos < pos_max))
    popt = np.array([np.nan] * 5)
    if len(i[0]):
        pos = pos[i]
        prof = prof[i]
        A0 = np.max(prof) - np.min(prof)
        offset0 = np.min(savgol_filter(prof, 21, 2))
        mu0 = pos[np.argmax(prof)]
        sigma0 = 1
        skew0 = 0
        if fit_func == 'Gauss4p':
            fitfunction = Gauss4p
            p0 = [A0, mu0, sigma0, offset0]
        elif fit_func == 'Gauss5p':
            fitfunction = Gauss5p
            p0 = [A0, mu0, sigma0, offset0, skew0]
        try:
            popt, pcov = curve_fit(fitfunction, pos, prof, p0=p0)
            if fit_func == 'Gauss4p': # added by Natalia, 14.10.2021
                errors = np.sqrt(np.array([pcov[0,0], pcov[1,1], pcov[2,2], pcov[3,3]])) 
            if fit_func == 'Gauss5p': # added from Natalia, 14.10.2021
                errors = np.sqrt(np.array([pcov[0,0], pcov[1,1], pcov[2,2], pcov[3,3], pcov[4,4]])) # Get the standard deviation of the parameters, square roors of the diagonal of the covariance
            if len(popt) == 4:
                popt = np.append(popt, [skew0])
        except (RuntimeError, ValueError):
            pass
    return popt, errors # return errors, added by Natalia, 14.10.2021

# From: https://gitlab.cern.ch/scripting-tools/injectors-beam-monitoring/-/blob/master/injectors_beam_monitoring/utils/analytical_functions.py
def Gauss3p(x, A, mu, sigma):
    return A * np.exp(-0.5 * (x - mu) ** 2 / sigma ** 2)


def Gauss4p(x, A, mu, sigma, offset):
    return Gauss3p(x, A, mu, sigma) + offset

def Gauss5p(x, A, mu, sigma, offset, skew):
    return Gauss4p(x, A, mu, sigma, offset) + skew * x

def getEmittance(sigma, beta_func, betagamma):
    return sigma ** 2 / beta_func * betagamma

def getEmittance_with_Dx(sigma, beta_func, betagamma, Dx, dpp):
	return (sigma**2-(dpp+Dx)**2)*betagamma/beta_func