"""
=========================================
Plotting light curves of rotating stars.
=========================================
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import kepler_data as kd
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel
from scipy.optimize import minimize


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

def neg_ln_like(p):
    gp.set_parameter_vector(p)
    return -gp.log_likelihood(y)

def grad_neg_ln_like(p):
    gp.set_parameter_vector(p)
    return -gp.grad_log_likelihood(y)


kepids = pd.read_csv("kepids.csv")
nstars = len(kepids.kepid.values)
filenames = ["/Users/ruthangus/.kplr/data/lightcurves/{0}"
                .format(str(kepids.kepid.values[i]).zfill(9))
                for i in range(nstars)]

# for i, filename in enumerate(filenames[:10]):
index = 8
filename = filenames[index]
x, y, yerr = kd.load_kepler_data(filename)
l, h = 2000, 5000
xfull, yfull, yerrfull = x[l:h]-x[l:h][0], y[l:h], yerr[l:h]
x, y, yerr = x[l:h:10]-x[l:h:10][0], y[l:h:10], yerr[l:h:10]

a, l, g, p = 1, 10, 10, 10
k = a * ExpSquaredKernel(l) * ExpSine2Kernel(g, p)
gp = george.GP(k)
gp.compute(x, yerr)

# optimize
result = minimize(neg_ln_like, gp.get_parameter_vector(),
                  jac=grad_neg_ln_like)
print(result)
print("A = ", result.x[0])
print("l = ", result.x[1])
print("gamma = ", result.x[2])
print("period = ", result.x[3])
gp.set_parameter_vector(result.x)
gp.compute(x, yerr)

x_pred = np.linspace(min(x), max(x), 500)
mu, var = gp.predict(y, x_pred, return_var=True)
mu *= 1e6
var *= 1e6

plt.figure(figsize=(20, 10))
plt.plot(xfull, yfull*1e6, "k.", ms=5)
plt.xlabel("$\mathrm{Time~[Days]}$")
plt.ylabel("$\mathrm{Flux~[parts~per~million]}$")
plt.xlim(0, max(x))
plt.savefig("lightcurves/{}_lc".format(kepids.kepid.values[index]))
plt.close()

color = "CornFlowerBlue"
plt.figure(figsize=(20, 10))
plt.plot(xfull, yfull*1e6, "k.", ms=5)
plt.plot(x_pred, mu, color=color, lw=2)
plt.fill_between(x_pred, mu - np.sqrt(var), mu + np.sqrt(var),
                 color=color, alpha=0.2)
plt.xlabel("$\mathrm{Time~[Days]}$")
plt.ylabel("$\mathrm{Flux~[parts~per~million]}$")
plt.xlim(0, max(x))
plt.savefig("lightcurves/{}_lc_gp".format(kepids.kepid.values[index]))
plt.close()
