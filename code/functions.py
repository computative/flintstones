from numpy import *
from matplotlib.pyplot import *

d = loadtxt("../resources/functions.txt")
fig = figure()
f, ax = subplots(2,sharex=True)
for n in range(6):
    x = d[n*500:500*(n+1)][:,0]
    y = d[n*500:500*(n+1)][:,1]
    z = d[n*500:500*(n+1)][:,2]
    ax[0].plot(x,z, label = ur"$n = %d$" % n)
    ax[1].plot(x,y, label = ur"$n = %d$" % n)
ax[0].set_ylim(-1,1)
ax[1].set_ylim(-50,50)
ax[1].set_xlabel(ur"$x \quad (\ \cdot \ )$")
ax[0].set_ylabel(ur"$\psi(x)\quad (\ \cdot\ )$")
ax[1].set_ylabel(ur"$H_n(x) \quad (\ \cdot\ )$")
ax[1].legend(fontsize=11)
ax[0].set_title(ur"HO basis elements and Hermite polynomials versus $x$")
savefig("../benchmark/functions.png")