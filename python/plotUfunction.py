# -*- coding: utf-8 -*-
"""
Simple script to plot the perpendicular integrals U(x,n) for n=0,1,2,3.
"""

import numpy as np
import matplotlib.pyplot as plt
import Hfunctions as hf

def plotU(index):
    names = [r"$U_0$",r"$\breve{U}_0$", r"$U_1$", r"$U_2$", r"$U_3$", r"$U_4$", r"$U_5$"]
    plt.figure()
    for harmonic in range(4):
        plt.plot(perpDistanceNormalized,UValues[index,harmonic,:],label=names[index]+', n='+str(harmonic),linewidth=3)
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$U(\xi)$')
    plt.legend()
    plt.xlim((0,upperLimit))
    plt.show()


resolution = 300
upperLimit = 3.0

perpDistanceNormalized = np.linspace(0,upperLimit,resolution)
UValues = np.zeros((7,4,resolution),dtype='float64')

for index in range(7):
    for harmonic in range(4):
        UValues[index,harmonic,:] = hf.getU(perpDistanceNormalized,index,harmonic)

plotU(0)
