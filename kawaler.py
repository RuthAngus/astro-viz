"""
Plot kawaler evolution.
"""

import numpy as np
import matplotlib.pyplot as plt

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def J2Omega(J, R, M):
    return J / (R**2 * M)


def dJdt(Omega, R, M, n=2.5, a=1, c=1):
    return Omega**(1 + 4*a*n/3) * R**(2-n) * M**(-n/3)


def dOdt(Omega, R, M):
    dJ_dt = dJdt(Omega, R, M)
    return J2Omega(dJ_dt, R, M)


def calc_periods(omega_init, R, M, delta_t, niter=100):
    """
    Evolve a star fwd in time and calculate its rotation period.
    """
    omega = omega_init*1
    ps, os = [np.zeros(niter) for i in range(2)]
    for i in range(niter):
        omega -= dOdt(omega, R, M) * delta_t
        os[i] = omega
        ps[i] = 2 * np.pi / omega
    return os, ps


if __name__ == "__main__":

    # periods = np.linspace(0, 20)
    # omegas = 2 * np.pi / periods
    # plt.plot(omegas, dOdt(omegas, 1, 1), ".")
    # plt.xlabel("omega")
    # plt.ylabel("dJ/dt")
    # plt.savefig("test")
    # plt.close()
    # assert 0

    R = np.array([1, 1, 1])
    M = np.array([1, 1, 1])
    period_init = np.array([2, 4, 4.5])
    omega_init = 2 * np.pi / period_init
    delta_t = 1e-2

    plt.figure(figsize=(16, 9))
    niter = 100
    for i in range(3):
        os, ps = calc_periods(omega_init[i], R[i], M[i], delta_t, niter=niter)
        # plt.plot(np.arange(niter) * delta_t * 1000, ps-.2*(period_init[i]))
        plt.plot(np.arange(niter) * delta_t * 500, ps)

    plt.xlabel("$\mathrm{Time~[Myr]}$")
    plt.ylabel("$P_{\mathrm{rot}}~\mathrm{[Days]}$")
    # plt.xlim(0, 500)
    plt.savefig("test")
    plt.close()
