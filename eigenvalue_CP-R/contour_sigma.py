# -*- coding: utf-8 -*-
"""
2025.3.5 created by kensuke ohtake

Eigenvalue analysis for CP-R
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime

# get datetime
DateTime = datetime.datetime.today().strftime("%Y%m%d%H%M%S")

# set parameters
m = 0.6 # mu
r = 1.0 # radius
v = 1.0 # adjust speed
lam = 1.0/(2.0*np.pi*r) # homogeneous manufactuirng population
ph = 1.0/(2.0*np.pi*r) # homogeneous agricultural population

def eigenv(X, Y, k):

    # X: tau mesh
    # Y: sigma mesh
    # k: frequency number (k)

    alp = (Y - 1.0) * X # alpha vector
    G = np.power(2.0*lam*((1.0-np.exp(-alp*np.pi*r))/alp), 1.0/(1.0-Y)) # homogeneous price index
    
    if k % 2 == 0:# when k is even number
        k = float(k)
        Z = (np.power(alp, 2.0)*np.power(r, 2.0)) / (np.power(k, 2.0) + np.power(alp, 2.0)*np.power(r, 2.0))
    else:# when k is odd number
        k = float(k)
        Z = (np.power(alp, 2.0)*np.power(r, 2.0)*(1.0 + np.exp(-alp*r*np.pi))) / ((np.power(k, 2.0) + np.power(alp, 2.0)*np.power(r, 2.0))*(1.0 - np.exp(-alp*r*np.pi)))

    Gamk = v * np.power(G, -m) * ((1.0-m*Z) * ((-(1.0/Y)*np.power(Z, 2.0)+(m/Y)*Z)/(1.0-(m/Y)*Z-((Y-1.0)/Y)*np.power(Z, 2.0))) + ((m*Z)/(Y-1.0)))
    return Gamk

tau_space = np.linspace(0.01, 2.5, 300)
sigma_space = np.linspace(1.01, 40.0, 1024)
ta, sig = np.meshgrid(tau_space, sigma_space)

frqs = [1,2,3,4,5,6]
labels = []
hs = []
fig, ax = plt.subplots()
i = 0
for k in frqs:
    Gamk = eigenv(ta, sig, k)
    cont = ax.contour(ta, sig, Gamk, [0], linewidths=1.5, colors=[matplotlib.cm.tab10(i)])
    lb = r'$k$ = {}'.format(k)
    labels.append(lb)
    h,_ = cont.legend_elements()
    hs.append(h[0])
    i += 1
    pass

ax.legend(hs, labels)
plt.gca().set_aspect('equal', adjustable='box')
#plt.grid()
ax.set_facecolor('0.95')
ax.set_aspect('auto', adjustable='box')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\sigma$')
plt.grid(linestyle=':')
plt.savefig('cpr_sigma_contours.png', format='png', dpi=300)
plt.show()
