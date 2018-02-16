"""
Plot stellar evolution tracks.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from isochrones.dartmouth import DartmouthModelGrid
import pandas as pd

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def zoom_out(df, llim):
    m = 10**df.log_L < llim + 1
    teff = 10**df.log_Teff.values[m]
    L = 10**df.log_L.values[m]
    star_age = df.star_age.values[m] * 1e-3
    R = 10**df.log_R.values[m]

    for i, age in enumerate(star_age[::5]):
        print(i, "of", len(star_age[::5]))
        m = age == star_age
        plt.clf()
        plt.figure(figsize=(16, 9))
        plt.title("Age = {0:.2f} Gyr".format(age*1e-6))
        plt.plot(teff, L, ".7", zorder=0)
        plt.scatter(teff[m], L[m], c=teff[m], vmin=3900, vmax=6100, s=R[m]*150,
                    zorder=1, cmap="plasma_r")
        plt.colorbar(label="$\mathrm{Effective~Temperature~[K]}$")
        plt.xlabel("$\mathrm{Effective~Temperature~[K]}$")
        plt.ylabel("$\mathrm{Luminosity}~[L_\odot]$")
        plt.xlim(6100, 3900)
        plt.ylim(.0, llim)
        plt.savefig("iso_movie/frame_{}".format(str(i).zfill(4)))
        plt.close()


def zoom_in(df, llim=2):
    m = (10**df.log_L < llim + 1) * (10**df.log_Teff < 6000) * \
        (5000 < 10**df.log_Teff)
    teff = 10**df.log_Teff.values[m]
    L = 10**df.log_L.values[m]
    star_age = df.star_age.values[m] * 1e-3
    R = 10**df.log_R.values[m]
    he = df.he_core_mass.values[m]

    for i, age in enumerate(star_age[::3]):
        print(i, "of", len(star_age[::3]))
        m = age == star_age
        plt.clf()
        plt.figure(figsize=(16, 9))
        plt.title("Age = {0:.2f} Gyr".format(age*1e-6))
        plt.plot(teff, L, ".7", zorder=0)
        plt.scatter(teff[m], L[m], c=teff[m], vmin=5500, vmax=6000, s=R[m]*150,
                    zorder=1, cmap="plasma_r")
        plt.colorbar(label="$\mathrm{Effective~Temperature~[K]}$")
        # plt.scatter(teff[m], L[m], c="w", vmin=3900, vmax=6100,
        #             s=(np.log(he[m])+5)*R[m]*300, zorder=2)
        plt.xlabel("$\mathrm{Effective~Temperature~[K]}$")
        plt.ylabel("$\mathrm{Luminosity}~[L_\odot]$")
        plt.xlim(6000, 5500)
        plt.ylim(.0, llim)
        plt.savefig("iso_movie_zoom_in/frame_{}".format(str(i).zfill(4)))
        plt.close()


def save_as_movie(movie_name, framerate=10, quality=25):
    """
    Make the movie file
    """
    os.system("/Applications/ffmpeg -r {0} -f image2 -s 1920x1080 -i "\
            "{1}/frame_%04d.png -vcodec libx264 -crf {2}  -pix_fmt "\
            "yuv420p {1}.mp4".format(framerate, movie_name, quality))


if __name__ == "__main__":
    path = "data/00100M.track.csv"
    df = pd.read_csv(path, skiprows=11)
    print(df.keys())
    # zoom_in(df)
    save_as_movie("iso_movie")
    save_as_movie("iso_movie_zoom_in")
