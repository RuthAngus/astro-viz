"""
==================
van Saders tracks
==================
Plotting van Saders model.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.constants as co
import astropy.units as u
import scipy.spatial as sps

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def logg(M, R):
    g = co.G * (M*co.M_sun) / (R*co.R_sun)**2
    g *= 1e2  # Convert to cgs
    return np.log10(g.value)


def period_age_figure(P, teff, age):
    plt.clf()
    plt.figure(figsize=(16, 9))
    plt.scatter(age, P, c=teff, s=20)
    plt.xlabel("$\mathrm{Age~[Gyr]}$")
    plt.ylabel("$\mathrm{Period~[Days]}$")
    plt.colorbar()
    plt.savefig("vs_tracks_period_age")


def period_teff_figure(P, teff, age):
    plt.clf()
    plt.figure(figsize=(16, 9))
    plt.scatter(teff, P, c=age, s=20)
    plt.xlabel("$T_\mathrm{eff}~[K]$")
    plt.ylabel("$\mathrm{Period~[Days]}$")
    plt.colorbar()
    plt.xlim(6250, 4000)
    plt.savefig("vs_tracks_period_teff")


def nearest(value, array):
    inds = (np.abs(value - array)).argmin()
    return inds, array(inds)


def animation(P, teff, age):
    t, p = 5770, 26
    print(np.shape(P), np.shape(teff), np.shape(age))
    for i, a in enumerate(np.linspace(0, 13.8, 10):
        inds, nearest_age = nearest(a, age)

        data = np.vstack((teff, P)).T
        tree = sps.cKDTree(data)
        dist, index = tree.query([t, p], 1)

        plt.clf()
        plt.figure(figsize=(16, 9))
        plt.scatter(teff, P, c=age, s=20, alpha=.1)
        plt.scatter(teff[inds], P[inds], c=age[inds], s=20)
        plt.xlabel("$T_\mathrm{eff}~[K]$")
        plt.ylabel("$\mathrm{Period~[Days]}$")
        plt.colorbar()#label="$\mathrm{Age~[Gyr]}$")
        plt.xlim(6250, 4000)
        plt.savefig("gyro_movie/frame_{}".format(str(j).zfill(4)))


if __name__ == "__main__":
    df = pd.read_csv("skeleton3_run_000.out", skiprows=172)
    print(df.keys())
    df["logg"] = logg(df.Mass_Msun.values, df.R_Rsun.values)
    m = (4000 < 10**df.log_Teff_K.values) & \
        (10**df.log_Teff_K.values < 6250) & (df.Age_Gyr > .1) & \
        (4.3 < df.logg.values)
    teff = 10**df.log_Teff_K.values[m]
    P = df.Prot_days.values[m]
    age = df.Age_Gyr.values[m]
    logg = df.logg.values[m]

    animation(P, teff, age)
