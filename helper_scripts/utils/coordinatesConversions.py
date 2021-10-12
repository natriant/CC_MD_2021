import numpy as np


def cmpt_normalised_coordinates(u, up, beta, alpha):
    # (u, up)--> (x,xp) or (y, yp), beta, alpha optic functions
    u_n = u / np.sqrt(beta)
    up_n = alpha*u/np.sqrt(beta) + np.sqrt(beta) * up
    return u_n, up_n


def cmpt_actions(u_n, up_n):
    J = (1 / 2) * (u_n ** 2 + up_n ** 2)
    return J
