# -*- coding: utf-8 -*-
"""
Simple script to plot the perpendicular integrals V(x,n) for n=0,1,2,3.
"""

import numpy as np
import matplotlib.pyplot as plt
import Hfunctions as hf

def plotV(index):
    names = ["V0","V0mod", "V1", "V2", "V3", "V4", "V5"]
    plt.figure()
    for harmonic in range(4):
        plt.plot(perpDistanceNormalized,VValues[index,harmonic,:],label=names[index]+', n='+str(harmonic),linewidth=3)
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$V(\xi)$')
    plt.legend()
    plt.xlim((0,upperLimit))
    plt.show()


resolution = 300
upperLimit = 3.0

perpDistanceNormalized = np.linspace(0,upperLimit,resolution)
VValues = np.zeros((7,4,resolution),dtype='float64')

for index in range(7):
    for harmonic in range(4):
        VValues[index,harmonic,:] = hf.getV(perpDistanceNormalized,index,harmonic)

plotV(0)
