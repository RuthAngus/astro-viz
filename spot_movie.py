'''
========================
Spot movie
========================
Make a movie of star spots.
'''

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)
# plt.rcParams['axes.facecolor'] = 'black'

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


def deg_to_rad(deg):
    return 2 * np.pi * deg / 360


def rad_to_deg(rad):
    return 360 * rad / (2*np.pi)


def convert_to_rads(lat_deg, lon_deg, lat_0_deg, lon_0_deg):
    # Convert to radians
    return deg_to_rad(lat_deg), deg_to_rad(lon_deg), \
        deg_to_rad(lat_0_deg), deg_to_rad(lon_0_deg)

def adjust_array_length(effective_lon):
    # Some fiddling needed because of finite pixel sizes.
    m = 0 < np.sin(effective_lon[i, :])
    pixel_array_len = len(effective_lon[i, :][m])
    if pixel_array_len < 89:
        j = 0
        while pixel_array_len < 89:
            # print(pixel_array_len, "pixel array too short!")
            m = (0 - j*.001) < np.sin(effective_lon[i, :])
            j += 1

    pixel_array_len = len(effective_lon[i, :][m])
    if pixel_array_len > 89:
        j = 0
        while pixel_array_len > 89:
            # print(pixel_array_len, "pixel array too long!")
            m = (0 + j*.001) < np.sin(effective_lon[i, :])
            pixel_array_len = len(effective_lon[i, :][m])
            j += 1
    return m


def integrated_light(lat_deg, lon_deg, Ti, lat_0_deg, lon_0_deg,
                     star_temp=4100):

    # Convert to radians
    lat, lon, lat_0, lon_0 = convert_to_rads(lat_deg, lon_deg, lat_0_deg,
                                             lon_0_deg)

    # Calculate effective latitude and longitude of each point on star.
    # This is lat/lon of each point + lat/lon of the star's orientation.
    effective_lon = lon + lon_0
    effective_lat = lat - lat_0

    # Calculate integrated flux
    flux, ws, Tis, max_flux, max_Ts = 0, 0, 0, 0, 0
    max_T = np.ones_like(Ti) * star_temp

    for i in range(np.shape(lon)[0]):  # For each row in lon:

        # Some fiddling needed because of finite pixel sizes.
        # m = adjust_array_length(effective_lon)

        # Select only Temperature pixels that are facing the observer.
        Tis += sum(Ti[i, :])

        # Select only weights that are facing the observer.
        # mw = np.sin(lon[i, :]) > 0
        # weights = np.sin(lat[i, :][mw]) * np.cos(lon[i, :][mw])
        # weights = np.sin(effective_lat[i, :]) * np.cos(effective_lon[i, :])
        # m = 0 < np.sin(lon[i, :])

        # weights = np.sin(lat[i, :]) * np.cos(lon[i, :])
        weights = np.sin(lat[i, :]) * np.cos(lon[i, :])
        m = 0 < np.sin(effective_lon[i, :]+np.pi/2)
        weights[m] = np.zeros(len(weights[m]))

        # Calculate the flux contribution for each row.
        # m = adjust_array_length(effective_lon)
        flux += sum(weights * Ti[i, :]) # [m])
        max_flux += sum(weights * max_T[i, :])  # [m])
        max_Ts += sum(max_T[i, :]) # [m])

    return flux, ws, Tis, max_flux


if __name__ == "__main__":

    nrows, ncols = 90, 180

    np.random.seed(123)
    # random data
    pts = 1 - 2 * np.random.rand(500, 3)
    l = np.sqrt(np.sum(pts**2, axis=1))  # Take the sum of the abs value.
    pts = pts / l[:, np.newaxis]

    relative_spot_flux = .001
    star_temp, spot_temp = 1./(nrows*ncols), relative_spot_flux/(nrows*ncols)
    T = star_temp * np.ones(500)  # np.random.rand(500)
    nspots = 10
    indices = np.random.choice(range(500), nspots)
    T[indices] = np.ones(nspots) * spot_temp

    # naive IDW-like interpolation on regular grid
    theta, phi, r = cart2sph(*pts.T)
    lon, lat = np.meshgrid(np.linspace(0, 360, ncols),
                           np.linspace(-90, 90, nrows))

    xg, yg, zg = sph2cart(lon, lat)
    Ti = np.zeros_like(lon)
    for r in range(nrows):
        for c in range(ncols):
            v = np.array([xg[r, c], yg[r, c], zg[r, c]])
            angs = np.arccos(np.dot(pts, v))
            idx = np.where(angs == 0)[0]
            if idx.any():
                Ti[r, c] = T[idx[0]]
            else:
                idw = 1 / angs**2 / sum(1 / angs**2)
                Ti[r, c] = np.sum(T * idw)

    nrotations = 2
    nframes = 50
    # nframes = 360*nrotations/20
    longitudes = np.linspace(360*nrotations, 0, nframes)
    fluxes, wgs, temps, max_flux = [], [], [], []
    for i, l in enumerate(longitudes):
        flux, weights, temps, max_flux = integrated_light(lat, lon, Ti, 10, l)
        fluxes.append(flux)
        wgs.append(weights)

    fluxes = np.array(fluxes) / max(fluxes)
    times = np.linspace(0, nrotations, len(longitudes))
    for i, l in enumerate(longitudes):
        plt.figure(figsize=(16, 9))
        plt.subplot(2, 1, 1)
        print(i, "of", len(longitudes))

        # set up map projection
        map = Basemap(projection='ortho', lat_0=10, lon_0=l)

        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0, 360, 30))
        map.drawparallels(np.arange(-90, 90, 30))

        # compute native map projection coordinates of lat/lon grid.
        x, y = map(lon, lat)

        # contour data over the map.
        # cs = map.contourf(x, y, Ti, 15, cmap="inferno_r")
        cs = map.contourf(lon, lat, Ti, 15, cmap="inferno", latlon=True)

        flux_ppm = fluxes  # *1e6
        plt.subplot(2, 1, 2)
        plt.plot(times[:i], flux_ppm[:i], ".")
        plt.xlim(0, 2)
        plt.ylim(min(flux_ppm)-1, max(flux_ppm)+1)
        plt.xlabel("$\\#~\mathrm{Full~Rotations}$")
        plt.ylabel("$\mathrm{Flux~[parts~per~million]}$")
        plt.savefig("spot_movie/frame_{}".format(str(i).zfill(4)))
        plt.close()
