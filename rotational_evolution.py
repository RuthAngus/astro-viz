"""
=======================================
Rotation evolution
=======================================
Making a video to show the evolution of stellar rotation periods over time as
a function of their mass.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import astropy.constants as co


def torque(R, M, B, period, Mdot=2e-14*1.989e30, K_1=1.3, K_2=.00506, m=.2177):
    """
    Calculate the torque acting on each star.
    Mdot is set to the Solar value.
    params:
    ------
    R: Radius (m)
    M: Mass (kg)
    B: Magnetic field strength (Gauss)
    Omega: Angular frequency = 2*pi/Period (1/seconds)
    """
    Omega = 2 * np.pi / period
    f = Omega*R**(3/2.) * (co.G*M)**-.5

    a = K_1 ** 2 / (2 * co.G)**m
    b = B**(4*m) * Mdot**(1 - 2*m)
    c = R**(5*m + 2) / M**m
    d = Omega / (K_2**2 + .5*f**2)**m
    return a * b * c * d


def spin_down_rate(tq, R, M, period, dt):
    """
    Calculate the new rotation period given a torque and a timestep.
    the rate of increasing rotation period for a given torque.
    torque = dL/dt
    dL = torque * dt
    dt: timestep (seconds)
    """
    dL = tq * dt

    Omega = 2*np.pi/period
    I = R**2 * M
    L = I * Omega
    new_L = L.value - dL.value
    new_period = I*2*np.pi/new_L
    return new_period


if __name__ == "__main__":
    # Make a plot of rotation period vs B-V colour. (later convert into mass)
    f_sun = 0.004  # Solar spin rate
    gamma_sun = 10**2.5  # (1e2 - 1e3)
    Mdot_sun = 2e-14  # Solar masses per year
    B_sun = 3*u.Gauss # Gauss (2-5)

    R = co.R_sun
    M = co.M_sun
    B = B_sun
    period = 26 * 24 * 60 * 60 * u.s
    Mdot = 2e-14 * 1.989e30 * u.kg/u.s
    t = torque(R, M, B, period, Mdot=Mdot, K_1=1.3, K_2=.00506, m=.2177)
    print(t.value)

    tq = t.value * u.m**2 * u.kg / u.s**2
    dt = 1e6 * 365.25 * 24 * 60 * 60 * u.s  # 1 million years
    new_period = spin_down_rate(t, R, M, period, dt)
    print(period.value/60/60/24, new_period.value/60/60/24)
