# Make some plots of exoplanet populations.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def period_radius(df):
    plt.clf()
    plt.plot(np.log10(df.pl_orbper), np.log10(df.pl_radj*radius_ratio), "k.",
            alpha=.5)
    plt.xlabel("$\log_{10}(\mathrm{Orbital~Period~[days]})$")
    plt.ylabel("$\mathrm{Planet~Radius~}[R_{\mathrm{Earth}}]$")
    # plt.xlim(0, 10)
    plt.subplots_adjust(left=.16, bottom=.15)
    plt.savefig("period_vs_radius")


def ra_dec(df):
    """
    Transform to xyz.
    """
    plt.rcParams['axes.facecolor'] = 'black'

    # Select only transiting planets.
    m = df.pl_kepflag.values == 1
    kepler = df.iloc[m]

    ra, dec = kepler.ra.values, kepler.dec.values
    c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    l = c_icrs.galactic.l.deg
    b = c_icrs.galactic.b.deg

    plt.clf()
    plt.figure(figsize=(20,20))
    plt.plot(l, b, "w.", ms=10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("l_vs_b")


def orbiting_planet():
    plt.rcParams['axes.facecolor'] = 'black'
    host_x, host_y = 0, 0

    a, period, slow_factor = 10, 10, 10
    random_phase = np.random.uniform(0, 2*np.pi)
    inc = 91

    nframes = 100
    x = a * np.cos(2*np.pi/(period*slow_factor) * np.arange(nframes)
                   + random_phase)
    z = a * np.sin(2*np.pi/(period*slow_factor) * np.arange(nframes)
                   + random_phase) * np.sin(inc)
    y = z
    # y = np.zeros_like(x)

    for i in range(nframes):
        plt.clf()
        plt.figure(figsize=(16, 9), dpi=300)
        if z[i] > 1:
            zdr_star, zdr_planet = 0, 1
        else:
            zdr_star, zdr_planet = 1, 0
        plt.plot(host_x, host_y, ".", color="w", ms=500, zorder=zdr_star)
        plt.plot(x[i], y[i], ".", color=".3", zorder=zdr_planet, ms=50)
        plt.xlim(-16, 16)
        plt.ylim(-9, 9)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig("orbit_movie/frame_{}".format(str(i).zfill(4)), dpi=300)
        plt.close()

    # Make the movie file
    import os
    # framerate = 30
    framerate = 20
    quality = 25
    os.system("/Applications/ffmpeg -r {0} -f image2 -s 1920x1080 -i "\
            "orbit_movie/frame_%04d.png -vcodec libx264 -crf {1}  -pix_fmt "\
            "yuv420p orbit_movie.mp4".format(framerate, quality))


if __name__ == "__main__":
    radius_ratio = 6.991e7/6.371e6
    df = pd.read_csv("planets.csv", skiprows=69)
    # ra_dec(df)

    orbiting_planet()
