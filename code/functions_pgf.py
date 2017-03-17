from numpy import *
# -*- coding: utf-8 -*-
from numpy import *
from scipy.stats import pareto, expon

import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
    "font.family": "serif",
    "font.serif": [],
    "font.size" : 11.0,
    "pgf.preamble": [
         r"\usepackage[utf8]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{cmbright}",
         r"\usepackage{newtxtext}",
         r"\usepackage{bm}",
         r"\usepackage{amsmath,amsthm}"
         ]
}
mpl.rcParams.update(pgf_with_pdflatex)
# blue, violet, green, brown, red ,orange
mpl.rcParams['axes.color_cycle'] = ['k', 'k','k', 'k', 'k', 'k']

import matplotlib.pyplot as plt
f, ax = plt.subplots(2,1,figsize=(2*0.8*4,2*0.8*3),sharex=True)

d = loadtxt("../resources/functions.txt")
#fig = figure()
#f, ax = subplots(2,sharex=True)
Style = [[1000],[1,1],[8,1,1,1,1,1,1,1],[4,4,2,2,1,1],[8,8],[8, 4, 2, 4, 2, 4]]
for n in range(6):
    x = d[n*500:500*(n+1)][:,0]
    y = d[n*500:500*(n+1)][:,1]
    z = d[n*500:500*(n+1)][:,2]
    ax[0].plot(x,z, label = ur"$n = %d$" % n,dashes=Style[n])
    ax[1].plot(x,y, label = ur"$n = %d$" % n,dashes=Style[n])
ax[0].set_ylim(-1,1)
ax[1].set_ylim(-50,50)
ax[1].set_xlabel(ur"$x \quad (\ \cdot \ )$")
ax[0].set_ylabel(ur"$\psi_n(x)\quad (\ \cdot\ )$")
ax[1].set_ylabel(ur"$H_n(x) \quad (\ \cdot\ )$")
ax[1].legend(fontsize=11)
ax[0].set_title(ur"Harmonic oscillator basis and Hermite polynomials output from class",fontsize=11)
plt.savefig("../benchmark/functions.pgf")
