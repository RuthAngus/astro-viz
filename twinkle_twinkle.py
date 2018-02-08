"""
=======================
Twinkling star movie
=======================
An animation of stars twinkling in the Kepler field
"""


import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import kepler_data as kd
import glob

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)
plt.rcParams['axes.facecolor'] = 'black'


def radec_to_lb(ra, dec):
    """
    Transform to ra and dec to l and b.
    """
    c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    l = c_icrs.galactic.l.deg
    b = c_icrs.galactic.b.deg
    return l, b


def get_ra_and_dec_from_headers(filenames):
    ra, dec = [], []
    for path in filenames:
        fname = glob.glob("{0}/*llc.fits".format(path))[0]
        with fits.open(fname) as hdul:
            hdr = hdul[0].header
        ra.append(hdr["RA_OBJ"])
        dec.append(hdr["DEC_OBJ"])
    return np.array(ra), np.array(dec)


def get_flux_array(nstars, ntimes, kepids, filenames):
    fluxes = np.zeros((nstars, ntimes))
    for i, star in enumerate(kepids.kepid.values[:nstars]):
        print(i, "of", len(kepids.kepid.values[:nstars]))
        _, flux, _ = kd.load_kepler_data(filenames[i])
        short_flux = flux[1000::10][:ntimes]

        # Scale the fluxes so they are between 0 and 1.
        short_flux = short_flux/(max(np.abs(short_flux))*2) + .5
        fluxes[i, :] = short_flux
    return fluxes


def plot_l_b(l, b, alphas, times):
    """
    l, b are 1d arrays containing the positions of the stars.
    Alphas are 2d arrays with the fluxes shape = ((nstars, times)).
    Each row is a different star and each column is a different epoch.
    t is the number of times or frames or fluxes.
    """
    for t in range(times):
        plt.clf()
        plt.figure(figsize=(20, 20))
        for i, star in enumerate(l):
            plt.plot(l[i], b[i], "w.", ms=10, alpha=alphas[i, t])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig("twinkle_movie/frame_{}".format(str(t).zfill(4)))


def save_as_movie(framerate=10, quality=25):
    """
    Make the movie file
    """
    os.system("/Applications/ffmpeg -r {0} -f image2 -s 1920x1080 -i "\
            "twinkle_movie/frame_%04d.png -vcodec libx264 -crf {1}  -pix_fmt "\
            "yuv420p twinkle_movie.mp4".format(framerate, quality))


if __name__ == "__main__":
    # df = pd.read_csv("planets.csv", skiprows=69)

    """
    Exoplanet host targets.
    """
    # m = df.pl_kepflag.values == 1
    # kepler = df.iloc[m]
    # ra, dec = kepler.ra.values, kepler.dec.values

    """
    Asteroseismology targets.
    """
    # kepids = pd.read_csv("kepids.csv")
    # nstars = len(kepids.kepid.values)
    # filenames = ["/Users/ruthangus/.kplr/data/lightcurves/{0}"
    #              .format(str(kepids.kepid.values[i]).zfill(9))
    #              for i in range(nstars)]

    # ntimes = 100
    # nstars = 525
    # ra, dec = get_ra_and_dec_from_headers(filenames)
    # fluxes = get_flux_array(nstars, ntimes, kepids, filenames)

    # l, b = radec_to_lb(ra, dec)
    # plot_l_b(l[:nstars], b[:nstars], fluxes, ntimes)
    save_as_movie(framerate=10, quality=25)
