# -*- coding: utf-8 -*-
"""
2025.3.2 created by Kensuke Ohtake

Eigenvalue analysis for QLLU-AD
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import datetime

# get datetime
DateTime = datetime.datetime.today().strftime("%Y%m%d%H%M%S")

# set parameters
m = 0.6 # mu
Lam = 1.0 # total manufacturing workers
Ph = 10.0 # total agricultural workers
r = 1.0 # radius
a = 0.5 # addvection coef.
d = 0.005 # diffusion coef.
lam = Lam/(2.0*np.pi*r) # homogeneous manufactuirng population
ph = Ph/(2.0*np.pi*r) # homogeneous agricultural population

def eigenv(X, Y, k):

    # X: tau mesh
    # Y: sigma mesh
    # k: frequency number (k)
    
    k = float(k)# 2025.4.10 updated

    alp = (Y - 1.0) * X # alpha vector
    
    if k % 2 == 0:# when k is even number
        Z = (np.power(alp, 2.0)*np.power(r, 2.0)) / (np.power(k, 2.0) + np.power(alp, 2.0)*np.power(r, 2.0))
    else:# when k is odd number
        Z = (np.power(alp, 2.0)*np.power(r, 2.0)*(1.0 + np.exp(-alp*r*np.pi))) / ((np.power(k, 2.0) + np.power(alp, 2.0)*np.power(r, 2.0))*(1.0 - np.exp(-alp*r*np.pi)))

    Gamk = (np.power(k, 2.0)/np.power(r, 2.0))*(a*m*Z*(-((ph+lam)/(Y*lam))*Z + ((2.0*Y-1.0)/(Y*(Y-1.0)))) - d)
    return Gamk

tau_space = np.linspace(0.01, 0.8, 256)
sigma_space = np.linspace(1.001, 7.0, 1025)
ta, sig = np.meshgrid(tau_space, sigma_space)

frqs = [1,2,3,4,5,6] # frequency numbers
for k in frqs:
    fig, ax = plt.subplots()
    Gamk = eigenv(ta, sig, k)
    norm = mcolors.TwoSlopeNorm(vmin=-0.15, vcenter=0, vmax=Gamk.max())
    levels = np.linspace(Gamk.min(), Gamk.max(), 256)
    cont = ax.contour(ta, sig, Gamk, [0], linewidths=1.5)
    contf = ax.contourf(ta, sig, Gamk, levels=levels, cmap='jet', norm=norm)
    #cmap='rainbow' 'bwr' 'coolwarm' 'seismic'    
    ticks = np.linspace(Gamk.min(), Gamk.max(), 9)
    fig.colorbar(contf, extend='max').set_ticks(ticks)
    
    ax.set_aspect('auto', adjustable='box')
    ax.set_title(r'k={}'.format(k))# 2025.4.17 updated

    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\sigma$')
    plt.savefig('sigma_heatmap_k_{}.png'.format(k), format='png', dpi=300)
    plt.show()
    pass
