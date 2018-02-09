"""
Make figures related to planet populations.
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plotpar = {'axes.labelsize': 30,
           'font.size': 25,
           'legend.fontsize': 25,
           'xtick.labelsize': 25,
           'ytick.labelsize': 25,
           'text.usetex': True}
plt.rcParams.update(plotpar)


df = pd.read_csv("planets.csv", skiprows=69)
m = df.pl_kepflag.values > 0
df = df.iloc[m]

# Fractional errorbars
perr1 = df.pl_orbpererr1.values / df.pl_orbper.values
perr2 = df.pl_orbpererr2.values / df.pl_orbper.values
rerr1 = df.pl_radjerr1.values / df.pl_radj.values
rerr2 = df.pl_radjerr2.values / df.pl_radj.values

# Log errorbars
logp = np.log10(df.pl_orbper.values)
logperr1 = logp*perr1
logperr2 = logp*perr2
logr = np.log10(df.pl_radj.values)
logrerr1 = logr*rerr1
logrerr2 = -logr*rerr2

m = (logperr1 < .5) * np.isfinite(logp) * np.isfinite(logr) * \
    np.isfinite(logperr1) * np.isfinite(logperr2) * np.isfinite(logrerr1) * \
    np.isfinite(logrerr2)

plt.clf()
plt.figure(figsize=(16, 9))
plt.errorbar(logp[m], logr[m], xerr=[logperr1[m], logperr2[m]],
             yerr=[logrerr1[m], logrerr2[m]], fmt="k.", ms=10, ecolor=".7")
plt.xlim(-1, 3.3)
plt.ylim(-1.7, .4)
plt.xlabel("$\log_{10}(\mathrm{Orbital~Period~[Days]})$")
plt.ylabel("$\log_{10}(\mathrm{Planet~Radius}[R_J])$")
plt.savefig("period_vs_radius")
