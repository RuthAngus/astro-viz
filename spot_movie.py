'''
========================
Spot movie
========================
Make a movie of star spots.
'''

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


def cart2sph(x, y, z):
    """
    Convert cartesian coords to spherical coords.
    """
    dxy = np.sqrt(x**2 + y**2)
    r = np.sqrt(dxy**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, dxy)
    theta, phi = np.rad2deg([theta, phi])
    return theta % 360, phi, r


def sph2cart(theta, phi, r=1):
    """
    Convert spherical coords to cartesian coords.
    """
    theta, phi = np.deg2rad([theta, phi])
    z = r * np.sin(phi)
    rcosphi = r * np.cos(phi)
    x = rcosphi * np.cos(theta)
    y = rcosphi * np.sin(theta)
    return x, y, z


if __name__ == "__main__":

    # random data
    pts = 1 - 2 * np.random.rand(500, 3)
    l = np.sqrt(np.sum(pts**2, axis=1))
    pts = pts / l[:, np.newaxis]
    T = 150 * np.ones(500)  # np.random.rand(500)
    nspots = 20
    indices = np.random.choice(range(500), nspots)
    T[indices] = np.ones(nspots) * 500

    # naive IDW-like interpolation on regular grid
    theta, phi, r = cart2sph(*pts.T)
    nrows, ncols = 90, 180
    lon, lat = np.meshgrid(np.linspace(0, 360, ncols),
                           np.linspace(-90, 90, nrows))

    xg, yg, zg = sph2cart(lon, lat)
    Ti = np.zeros_like(lon)
    for r in range(nrows):
        for c in range(ncols):
            v = np.array([xg[r, c], yg[r, c], zg[r, c]])
            angs = np.arccos(np.dot(pts, v))
            idx = np.where(angs == 0)[0]
            # if idx.any():
            #     Ti[r, c] = T[idx[0]]
            # else:
            idw = 1 / angs**2 / sum(1 / angs**2)
            Ti[r, c] = np.sum(T * idw)

    nframes = 360/10
    longitudes = np.linspace(360, 0, nframes)

    for i, l in enumerate(longitudes):
        print(i, "of", len(longitudes))
        # set up map projection
        map = Basemap(projection='ortho', lat_0=10, lon_0=l)
        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0, 360, 30))
        map.drawparallels(np.arange(-90, 90, 30))
        # compute native map projection coordinates of lat/lon grid.
        x, y = map(lon, lat)
        # contour data over the map.
        cs = map.contourf(x, y, Ti, 15, cmap="inferno_r")
        plt.savefig("spot_movie/frame_{}".format(str(i).zfill(4)))
        plt.close()
