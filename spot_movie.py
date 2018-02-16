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

    effective_lon = deg_to_rad(lon_deg + lon_0_deg)
    effective_lat = deg_to_rad(lat_deg + lat_0_deg)

    # Convert to radians
    lat, lon, lat_0, lon_0 = convert_to_rads(lat_deg, lon_deg, lat_0_deg,
                                             lon_0_deg)

    # Calculate effective latitude and longitude of each point on star.
    # This is lat/lon of each point + lat/lon of the star's orientation.

    # Calculate integrated flux
    flux, ws, Tis, max_flux, max_Ts = 0, 0, 0, 0, 0
    max_T = np.ones_like(Ti) * star_temp
    weight_array = np.zeros_like(Ti)

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

        # weights = np.cos(effective_lat[i, :]) * np.cos(effective_lon[i, :])
        weights = np.cos(lat[i, :]) * np.cos(lon[i, :])
        m = (0 < np.cos(effective_lon[i, :])) #* \
            # (0 < np.cos(effective_lat[i, :]))
        weights[m] = np.zeros(len(weights[m]))
        weights = np.abs(weights)
        weight_array[i, :] = weights

        # Calculate the flux contribution for each row.
        # m = adjust_array_length(effective_lon)
        flux += sum(weights * Ti[i, :]) # [m])
        max_flux += sum(weights * max_T[i, :])  # [m])
        max_Ts += sum(max_T[i, :]) # [m])

#         plt.figure(figsize=(16, 9))
#         plt.plot(effective_lon[i, :]/np.pi, weights)
#         plt.xlabel("longitude [pi]")
#         plt.ylabel("weight")
#         plt.savefig("spot_movie/lon_{}".format(str(i).zfill(4)))
#         plt.close()

        # plt.figure(figsize=(16, 9))
        # plt.plot(effective_lon[i, :], Ti[i, :])
        # plt.ylabel("temp")
        # plt.savefig("spot_movie/temp_{}".format(str(i).zfill(4)))
        # plt.close()

        # plt.plot(lon[i, :], weights)
        # plt.ylabel("weights")
        # plt.ylim(-1, 1)
        # plt.savefig("spot_movie/weights_{}".format(str(i).zfill(4)))
        # plt.close()

    # plt.figure(figsize=(16, 9))
    # plt.plot(effective_lat[:, 45]/np.pi, weight_array[:, 45])
    # plt.xlabel("latitude [pi]")
    # plt.ylabel("weight")
    # # plt.savefig("spot_movie/lat_{}".format(str(i).zfill(4)))
    # plt.savefig("lats")
    # plt.close()
    # assert 0

    return flux, ws, Tis, max_flux, weight_array


def add_spot(x, y, z, star_radius, spot_phi, spot_theta, spot_radius):
    spot_center_uv = np.array([spot_phi,spot_theta])
    spot_center_xyz = np.array([np.outer(np.cos(spot_center_uv[0]),np.sin(spot_center_uv[1])),
                                np.outer(np.sin(spot_center_uv[0]),
                                         np.sin(spot_center_uv[1])),
                                np.outer(1,np.cos(spot_center_uv[1]))])
    x_sep = x/star_radius - spot_center_xyz[0]
    y_sep = y/star_radius - spot_center_xyz[1]
    z_sep = z/star_radius - spot_center_xyz[2]
    chord_length = np.sqrt(x_sep**2 + y_sep**2 + z_sep**2)
    central_angle = 2 * np.arcsin(chord_length/2)
    great_circle_distance = star_radius * central_angle
    spot = great_circle_distance < spot_radius
    return spot


if __name__ == "__main__":

    # nrows, ncols = 180, 360
    nrows, ncols = 90, 180

    np.random.seed(123)
    # random data
    pts = 1 - 2 * np.random.rand(500, 3)
    l = np.sqrt(np.sum(pts**2, axis=1))  # Take the sum of the abs value.
    pts = pts / l[:, np.newaxis]

    relative_spot_flux = 0
    star_temp, spot_temp = 1e8/(nrows*ncols), relative_spot_flux#/(nrows*ncols)
    T = star_temp * np.ones(500)  # np.random.rand(500)
    nspots = 1
    indices = np.random.choice(range(500), nspots)
    # T[indices] = np.ones(nspots) * spot_temp

    # naive IDW-like interpolation on regular grid
    theta, phi, r = cart2sph(*pts.T)
    lon, lat = np.meshgrid(np.linspace(0, 360, ncols),
                           np.linspace(-90, 90, nrows))
                           # np.linspace(0, 180, nrows))

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

    star_radius, spot_theta, spot_phi, spot_radius = 1, 0, np.pi/2, .5
    spot = add_spot(xg, yg, zg, star_radius, spot_theta, spot_phi, spot_radius)
    for i, row in enumerate(Ti[:, 0]):
        m = (0 < lon[i, :]) * (lon[i, :] < 30)
        Ti[i, :][spot[i, :]] = np.zeros(len(Ti[i, :][spot[i, :]]))

    plt.imshow(Ti)
    plt.colorbar()
    plt.savefig("temp_map")
    plt.close()

    nrotations = 2
    nframes = 50
    # nframes = 360*nrotations/20
    longitudes = np.linspace(360*nrotations, 0, nframes)
    fluxes, wgs, temps, max_flux = [], [], [], []
    for i, l in enumerate(longitudes):
        flux, weights, temps, max_flux, weight_array = \
            integrated_light(lat, lon, Ti, 0, l)
        fluxes.append(flux)
        wgs.append(weights)

    plt.imshow(weight_array)
    plt.colorbar()
    plt.savefig("weight_map")
    plt.close()

    # fluxes = np.array(fluxes) / max(fluxes)
    fluxes = np.array(fluxes) / np.var(fluxes)
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
        plt.plot(times[:i], -flux_ppm[:i], ".")
        plt.xlim(0, 2)
        # plt.ylim(min(flux_ppm)-1, max(flux_ppm)+1)
        plt.xlabel("$\\#~\mathrm{Full~Rotations}$")
        plt.ylabel("$\mathrm{Flux~[parts~per~million]}$")

        plt.savefig("spot_movie/frame_{}".format(str(i).zfill(4)))
        plt.close()
