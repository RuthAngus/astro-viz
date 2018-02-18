"""
plot the habitable zones of the Sun and an M dwarf.
"""


import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as co
import pandas as pd

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)
plt.rcParams['axes.facecolor'] = 'black'


def equilibrium_temperature(T_star, R_star, albedo, semi_major_axis):
    """
    Calculate the equilibrium temperature of a planet.
    """
    return T_star * (1 - albedo)**.25 * (R_star/(2*semi_major_axis))**.5


def HZ(T_star, R_star, albedo, T_p):
    """
    Calculate the Hz
    """
    return .5 * R_star * T_star**2 * (1 - albedo)**.5 / T_p**2


def plot(inner, outer, teff, r, a, min_teff, max_teff, logage, name, fnumber):
    y = np.linspace(-.05, .05, 100)
    x_inner = (inner**2 - y**2)**.5
    x_outer = (outer**2 - y**2)**.5

    plt.figure(figsize=(16, 4), dpi=300)
    plt.scatter(0, 0, c=min_teff, s=r, vmin=min_teff, vmax=max_teff,
                cmap="YlOrRd_r")
    plt.scatter(0, 0, c=teff, s=r, vmin=min_teff, vmax=max_teff,
                cmap="YlOrRd_r")
    plt.plot(a, 0, "o", color="mediumaquamarine", ms=10)
    plt.fill_betweenx(y, x_inner, x_outer, alpha=.5)
    # plt.fill_betweenx(y, 0, x_inner, alpha=.5, color="red")
    # plt.fill_betweenx(y, x_outer, 10, alpha=.5, color="blue")
    plt.title("Age = {0:.3f} Gyr".format(np.exp(logage)))
    plt.colorbar(label="$\mathrm{Effective~Temperature~[K]}$")
    plt.ylim(-.05, .05)
    if a < 1:
        plt.ylim(-.015, .015)
    plt.xlim(0, 2*a)
    plt.savefig("{}/frame_{}".format(name, str(fnumber).zfill(4)), dpi=300)
    plt.close()


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def HZ_movie(df, min_teff, max_teff, max_L, a, name, interval=10):
    # Select a small section of the HR diagram.
    m = (10**df.log_L < max_L) * (10**df.log_Teff < max_teff) * \
        (min_teff < 10**df.log_Teff) * (df.star_age*1e-9 < 13.8)
    teff = 10**df.log_Teff.values[m]
    L = 10**df.log_L.values[m]
    star_age = df.star_age.values[m] * 1e-9
    R = 10**df.log_R.values[m]
    he = df.he_core_mass.values[m]

    # for i, logage in enumerate(np.log(star_age[::interval])):
    nsteps = 50
    logages = np.log(np.linspace(0, 13.8, nsteps))
    for i, la in enumerate(logages):
        print(i, "of", nsteps)
        logage = find_nearest(np.log(star_age), la)
        star = logage == np.log(star_age)
        outer_HZ = HZ(teff[star], R[star]*co.R_sun, a, T_outer)
        inner_HZ = HZ(teff[star], R[star]*co.R_sun, a, T_inner)
        inner = inner_HZ/co.au
        outer = outer_HZ/co.au
        print(inner)
        plot(inner[0], outer[0], teff[star][0], 4000*R[star][0], a, min_teff,
             max_teff, logage, name, i)


if __name__ == "__main__":
    T_sun = 5777
    R_sun = co.R_sun
    a_earth = .3
    T_inner, T_outer = 373, 250
    print(equilibrium_temperature(T_sun, R_sun, a_earth, co.au))

    min_teff, max_teff, max_L = 3000, 7000, 50
    # Load isochrones
    path = "data/00100M.track.csv"
    df = pd.read_csv(path, skiprows=11)
    # HZ_movie(df, min_teff, max_teff, max_L, a_earth, "sun_HZ")

    min_teff, max_teff, max_L, a = 2000, 5000, 10, .025
    path = "data/00010M.track.csv"
    df = pd.read_csv(path, skiprows=11)
    HZ_movie(df, min_teff, max_teff, max_L, a, "mdwarf_HZ", interval=3)
