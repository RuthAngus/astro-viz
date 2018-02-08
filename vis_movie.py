# Make frames of exoplanet orbits.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)
plt.rcParams['axes.facecolor'] = 'black'


def get_orbit(a, radius, period, frames, xcenter=0, ycenter=0,
              slow_factor=30):
    """
    nframes = snapshots per orbit.
    if period long snapshots per orbit goes up.
    if period short, snapshots per orbit goes down.
    snapshots per orbit (nframes) is proportional to the period and
    1/angular frequency.
    """
    random_phase = np.random.uniform(0, 2*np.pi)
    x = a * np.cos(np.pi/(period*slow_factor) * np.arange(nframes)
                   + random_phase)
    y = a * np.sin(np.pi/(period*slow_factor) * np.arange(nframes)
                   + random_phase)
    return x, y

if __name__ == "__main__":

    # Read in the data
    df = pd.read_csv("planets.csv", skiprows=69)

    # Select only transiting planets.
    m = []
    for i, method in enumerate(list(df.pl_discmethod.values)):
        if method == "Transit":
            m.append(i)
    df = df.iloc[np.array(m)]

    R_j = 6.991e7
    radius_ratio = R_j/6.371e6  # Convert to Earth radii
    nframes = 800
    assert nframes % 2 == 0, "nframes must be an even number"
    nstar = 800

    a_s = df.pl_orbsmax.values * 1.496e11 / R_j
    radii = df.pl_radj.values * radius_ratio
    periods = df.pl_orbper.values
    xpositions = np.zeros((nstar, nframes))
    ypositions = np.zeros((nstar, nframes))
    for star in range(nstar):
        x, y = get_orbit(a_s[star], radii[star], periods[star], nframes)
        xpositions[star, :] = x
        ypositions[star, :] = y

    lim_min = -np.logspace(-2, np.log10(.5), nframes)
    lim_max = np.logspace(-2, np.log10(.5), nframes)
    for pos in range(nframes):
        plt.clf()
        plt.figure(figsize=(16, 9))
        for star in range(nstar):
            plt.scatter(xpositions[star, pos], ypositions[star, pos],
                        s=radii[star]*1./(lim_max[pos]-lim_min[pos]),
                        c="w", alpha=.5)
        plt.xlim(lim_min[pos], lim_max[pos])
        plt.ylim(lim_min[pos], lim_max[pos])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig("planet_cloud_movie/frame_{}".format(str(pos).zfill(4)))
        plt.close()

    # Make the movie file
    import os
    # framerate = 30
    framerate = 20
    quality = 25
    os.system("/Applications/ffmpeg -r {0} -f image2 -s 1920x1080 -i "\
            "planet_cloud_movie/frame_%04d.png -vcodec libx264 -crf {1} "\
             "-pix_fmt yuv420p planet_cloud_movie.mp4".format(framerate,
                                                              quality))
