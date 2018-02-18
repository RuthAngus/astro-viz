"""
plot the habitable zones of the Sun and an M dwarf.
"""


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.facecolor'] = 'black'
import astropy.constants as co


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


if __name__ == "__main__":
    T_sun = 5777
    R_sun = co.R_sun
    a_earth = .3
    T_inner, T_outer = 373, 250

    outer_HZ = HZ(T_sun, R_sun, a_earth, T_outer)
    inner_HZ = HZ(T_sun, R_sun, a_earth, T_inner)
    print(inner_HZ/co.au, outer_HZ/co.au)
    inner = inner_HZ/co.au
    outer = outer_HZ/co.au

    print(equilibrium_temperature(T_sun, R_sun, a_earth, co.au))

    y = np.linspace(-.05, .05, 100)
    x_inner = (inner**2 - y**2)**.5
    x_outer = (outer**2 - y**2)**.5

    # Plot the Sun, earth and habitable zone.
    fig = plt.figure(figsize=(16, 4), dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(0, 0, "o", color="Gold", ms=100)
    ax.plot(1, 0, "o", color="mediumaquamarine", ms=10)
    ax.fill_betweenx(y, x_inner, x_outer, alpha=.5)
    ax.set_ylim(-.05, .05)
    ax.set_xlim(0, 2)
    plt.savefig("sun_HZ", dpi=300)
    plt.close()
