"""
=======================
Twinkling star movie
=======================
An animation of stars twinkling in the Kepler field
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import kepler_data as kd

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)
plt.rcParams['axes.facecolor'] = 'black'


def ra_dec(ra, dec):
    """
    Transform to ra and dec to l and b.
    """

    # Select only transiting planets.
    m = df.pl_kepflag.values == 1
    kepler = df.iloc[m]

    ra, dec = kepler.ra.values, kepler.dec.values
    c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    l = c_icrs.galactic.l.deg
    b = c_icrs.galactic.b.deg
    return l, b


def plot_l_b(l, b, alphas, times):
    """
    l, b and alphas are 2d arrays with shape ((nstars, times)).
    Each row is a different star and each column is a different epoch.
    """
    for t in times:
        plt.clf()
        plt.figure(figsize=(20,20))
        for i, star in enumerate(l):
            plt.plot(l[i, t], b[i, t], "w.", ms=10, alpha=alphas[i, t])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig("twinkle/frame_{}".format(t))


if __name__ == "__main__":
    # df = pd.read_csv("planets.csv", skiprows=69)
    # ra_dec(df)

    kepids = pd.read_csv("kepids.csv")
    nstars = len(kepids.kepid.values)
    filenames = ["/Users/ruthangus/.kplr/data/lightcurves/{0}"
                 .format(str(kepids.kepid.values[i]).zfill(9))
                 for i in range(nstars)]

    ntimes = 100
    fluxes = np.zeros((nstars, ntimes))
    for i, star in enumerate(kepids.kepid.values):
        _, flux, _ = kd.load_kepler_data(filenames[0])
        fluxes[i, :] = flux[::100][:100]
