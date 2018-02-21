'''
========================
Spot movie
========================
Make a movie of star spots.
'''

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from matplotlib.colors import LightSource
from matplotlib import cm

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)
plt.rcParams['axes.facecolor'] = 'black'


def deg_to_rad(deg):
    return 2 * np.pi * deg / 360


def rad_to_deg(rad):
    return 360 * rad / (2*np.pi)


def convert_to_rads(lat_deg, lon_deg, lat_0_deg, lon_0_deg):
    # Convert to radians
    return deg_to_rad(lat_deg), deg_to_rad(lon_deg), \
        deg_to_rad(lat_0_deg), deg_to_rad(lon_0_deg)


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

        # weights = np.cos(effective_lat[i, :]) * np.cos(effective_lon[i, :])
        weights = np.sin(lat[i, :]) * np.cos(lon[i, :])
        m = (0 < np.cos(effective_lon[i, :])) #* \
            # (0 < np.cos(effective_lat[i, :]))
        weights[m] = np.zeros(len(weights[m]))
        weights = np.abs(weights)
        weight_array[i, :] = weights

        # Calculate the flux contribution for each row.
        flux += sum(weights * Ti[i, :]) # [m])
        max_flux += sum(weights * max_T[i, :])  # [m])
        max_Ts += sum(max_T[i, :]) # [m])

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


def convert_to_rads(lat_deg, lon_deg, lat_0_deg, lon_0_deg):
    # Convert to radians
    return deg_to_rad(lat_deg), deg_to_rad(lon_deg), \
        deg_to_rad(lat_0_deg), deg_to_rad(lon_0_deg)


def integrated_light(lat_deg, lon_deg, Ti, lat_0_deg, lon_0_deg,
                     star_temp=4100):

    # Calculate effective latitude and longitude of each point on star.
    # This is lat/lon of each point + lat/lon of the star's orientation.
    effective_lon = deg_to_rad(lon_deg + lon_0_deg)
    effective_lat = deg_to_rad(lat_deg + lat_0_deg)

    # Convert to radians
    lat, lon, lat_0, lon_0 = convert_to_rads(lat_deg, lon_deg, lat_0_deg,
                                             lon_0_deg)

    # Calculate integrated flux
    flux, ws, Tis, max_flux, max_Ts = 0, 0, 0, 0, 0
    max_T = np.ones_like(Ti) * star_temp
    weight_array = np.zeros_like(Ti)

    for i in range(np.shape(lon)[0]):  # For each row in lon:
        # Select only Temperature pixels that are facing the observer.
        Tis += sum(Ti[i, :])

        weights = np.sin(effective_lat[i, :]) * np.cos(effective_lon[i, :])
        m = (0 < np.cos(effective_lon[i, :]))
        weights[m] = np.zeros(len(weights[m]))
        weights = np.abs(weights)
        weight_array[i, :] = weights

        # Calculate the flux contribution for each row.
        flux += sum(weights * Ti[i, :])
        max_flux += sum(weights * max_T[i, :])
        max_Ts += sum(max_T[i, :])

    return flux, ws, Tis, max_flux, weight_array


if __name__ == "__main__":

    np.random.seed(1234)
    star_radius = 1

    phi = np.linspace(0, 2*np.pi, 512*2)
    theta = np.linspace(0, np.pi, 256*2)

    r = star_radius

    x = r*np.outer(np.cos(phi), np.sin(theta))
    y = r*np.outer(np.sin(phi), np.sin(theta))
    z = r*np.outer(np.ones(np.size(phi)), np.cos(theta))

    nspots = 10
    spot_theta = np.random.uniform(0, 2*np.pi, nspots)
    spot_phi = np.random.uniform(-np.pi/2, np.pi/2, nspots)
    spot_radius = np.exp(np.random.uniform(-4, -2, nspots))

    # Add spots
    c = np.empty(z.shape , dtype="str")
    Ti = np.empty(z.shape)
    Ti[:] = 1
    c[:] = "w"
    for s in range(nspots):
        spot = add_spot(x, y, z, r, spot_phi[s], spot_theta[s],
                        spot_radius[s])
        c[spot] = "k"
        Ti[spot] = 0

    # Calculate light curve.
    nrotations = 3
    nframes = 100
    longitudes = np.linspace(360*nrotations, 0, nframes)
    fluxes, wgs, temps, max_flux = [], [], [], []
    lon, lat, _ = cart2sph(x, y, z)
    for i, l in enumerate(longitudes):
        flux, weights, temps, max_flux, weight_array = \
            integrated_light(lat, lon, Ti, 0, l)
        fluxes.append(flux)
        wgs.append(weights)
    fluxes = np.array(fluxes) / np.var(fluxes)
    times = np.linspace(0, nrotations, len(longitudes))

    for i, azimuth in enumerate(longitudes):
        print(i)

        # print(type(light))
        # illuminated_surface = light.shade(z.shape)
        # illuminated_surface = light.shade(z.shape, cmap=cm.coolwarm)

        light = LightSource(90, azimuth)
        c = np.empty(z.shape , dtype="str")
        c[:] = "w"
        cz = np.random.randn(np.shape(z)[0], np.shape(z)[1]) + 4100

        # Add spots
        Ti = np.empty(z.shape)
        Ti[:] = 1
        for s in range(nspots):
            spot = add_spot(x, y, z, r, spot_phi[s], spot_theta[s],
                            spot_radius[s])
            c[spot] = "k"
            Ti[spot] = 0
            cz[spot] = 2000

        c = light.shade(cz, cmap=cm.inferno)

        fig = plt.figure(figsize=(16, 9), dpi=300)

        gs = gridspec.GridSpec(4, 3)
        # ax = fig.add_subplot(211, projection='3d', aspect='equal')
        ax1 = fig.add_subplot(gs[:3, :], projection='3d', aspect='equal')
        ax1.set_axis_off()
        ax1.view_init(elev=0, azim=azimuth)

        ax1.plot_surface(x, y, z, facecolors=c, rstride=1, cstride=1)

        # ax = fig.add_subplot(212)
        ax2 = fig.add_subplot(gs[-1:, :])
        flux_ppm = fluxes  # *1e6
        ax2.plot(times[:i], -flux_ppm[:i], "w.", ms=10)
        ax2.set_xlim(0, nrotations)
        ax2.set_ylim(min(-flux_ppm)-.01, max(-flux_ppm)+.01)
        ax2.set_xlabel("$\\#~\mathrm{Full~Rotations}$")
        # ax2.set_ylabel("$\mathrm{Flux~[parts~per~million]}$")
        ax2.set_ylabel("$\mathrm{Flux}$")

        plt.savefig("spot_movie/frame_{}".format(str(i).zfill(4)), dpi=300)
