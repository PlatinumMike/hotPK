# -*- coding: utf-8 -*-
"""
Simple script to plot the perpendicular integrals H(x,n) for n=0,1,2,3
"""

import numpy as np
import matplotlib.pyplot as plt
import Hfunctions as hf

def plotH(index,jacobian=False):
    names = [r"$H_0$",r"$\breve{H}_0$", r"$H_1$", r"$H_2$", r"$H_3$", r"$H_4$", r"$H_5$"]
    plt.figure()
    for harmonic in range(4):
        plt.plot(perpDistanceNormalized,HValues[index,harmonic,:],label=names[index]+', n='+str(harmonic),linewidth=3)
    plt.xlabel(r'$\xi$')
    if jacobian:
        plt.ylabel(r'$\xi H(\xi)$')
    else:
        plt.ylabel(r'$H(\xi)$')
    plt.legend()
    plt.xlim((0,upperLimit))
    plt.show()

resolution = 300
upperLimit = 3.0
jacobian = True

perpDistanceNormalized = np.linspace(0,upperLimit,resolution)
HValues = np.zeros((7,4,resolution),dtype='float64')

#evaluate on some grid
for index in range(7):
    for harmonic in range(4):
        HValues[index,harmonic,:] = hf.getH(perpDistanceNormalized,index,harmonic,jacobian)

plotH(1,jacobian)
plt.tight_layout()