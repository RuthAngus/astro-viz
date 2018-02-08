"""
=========================================
Plotting light curves of rotating stars.
=========================================
"""

import numpy
import matplotlib.pyplot as plt
import pandas as pd
import kepler_data as kd

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)

"""
Asteroseismology targets.
"""

kepids = pd.read_csv("kepids.csv")
nstars = len(kepids.kepid.values)
filenames = ["/Users/ruthangus/.kplr/data/lightcurves/{0}"
                .format(str(kepids.kepid.values[i]).zfill(9))
                for i in range(nstars)]

# for i, filename in enumerate(filenames[:10]):
index = 8
filename = filenames[index]
x, y, yerr = kd.load_kepler_data(filename)
plt.clf()
plt.figure(figsize=(20, 10))
l, h = 2000, 5000
x, y = x[l:h]-x[l:h][0], y[l:h]
plt.plot(x, y*1e6, "k.", ms=5)
plt.xlabel("$\mathrm{Time~[Days]}$")
plt.ylabel("$\mathrm{Flux~[parts~per~million]}$")
plt.xlim(0, max(x))
plt.savefig("lightcurves/{}_lc".format(kepids.kepid.values[index]))
plt.close()
