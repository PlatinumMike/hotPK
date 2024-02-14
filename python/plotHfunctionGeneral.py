# -*- coding: utf-8 -*-
"""
Simple script to plot the perpendicular integrals H(x,n) for integer n
"""

import numpy as np
import matplotlib.pyplot as plt
import HfunctionGeneral as hfg


def plotH(index, jacobian=False):
    names = ["H0", "H0mod", "H1", "H2", "H3", "H4", "H5"]
    plt.figure()
    for harmonic in range(Nharmonics):
        plt.plot(
            perpDistanceNormalized,
            HValues[index, harmonic, :],
            label=names[index] + ", n=" + str(harmonic),
            linewidth=3,
        )
    plt.xlabel(r"$\xi$")
    if jacobian:
        plt.ylabel(r"$\xi H(\xi)$")
    else:
        plt.ylabel(r"$H(\xi)$")
    plt.legend()
    plt.xlim((0, upperLimit))
    plt.show()


resolution = 300
upperLimit = 3.0
jacobian = True
Nharmonics = 8

perpDistanceNormalized = np.linspace(0, upperLimit, resolution)
HValues = np.zeros((7, Nharmonics, resolution), dtype="float64")

# evaluate on some grid
for index in range(7):
    for harmonic in range(Nharmonics):
        for point in range(resolution):
            HValues[index, harmonic, point] = float(
                hfg.getH(perpDistanceNormalized[point], index, harmonic, jacobian)
            )

plotH(1, jacobian)
