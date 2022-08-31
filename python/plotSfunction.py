# -*- coding: utf-8 -*-
"""
Simple script to plot the parallel integrals S(mu).
"""

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
import Sfunctions as sf

mp.mp.dps=100
plt.rcParams.update({'font.size': 16})

def castToNumpy(mpArray):
    return np.array(mpArray.tolist(),dtype='complex128')[0,:]

def cutZero(arr):
    length = arr.size
    if length%2==1 and length>2:
        #odd count, so assuming 0 is in the middle
        return np.concatenate((arr[:length//2],arr[length//2+1:]))
    else:
        print("Error: length is even")


resolution = 301 #odd so that we have a point exactly at mu=0.
eps = 1 #resonance charge sign: -1 or 1.

mu = mp.linspace(-3,3,resolution)

S1 = mp.matrix(1,resolution)
S2 = mp.matrix(1,resolution)
S3 = mp.matrix(1,resolution)

for i in range(resolution):
    S1[i] = sf.getS1(eps, mu[i])
    S2[i] = sf.getS2(eps, mu[i])
    S3[i] = sf.getS3(eps, mu[i])

#convert for plotting. Sure, mpmath can plot natively as well, but I'd like to export the values as well, so best to convert to double precision.
munumpy = np.array(mu,dtype='float64')
S1numpy = castToNumpy(S1)
S2numpy = castToNumpy(S2)
S3numpy = castToNumpy(S3)

plt.figure()
plt.plot(munumpy,S1numpy.real,label=r'$S_1$',linewidth=2)
plt.plot(munumpy,S2numpy.real,label=r'$S_2$',linewidth=2)
plt.plot(munumpy,S3numpy.real,label=r'$S_3$',linewidth=2)
plt.xlim(-3,3)
plt.ylim(-6.5,6.5)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\mathrm{Re}(S)$')
plt.legend()
plt.tight_layout()

#Im part has some trouble at mu=0, so eliminate one point
muCut = cutZero(munumpy)
S1cut = cutZero(S1numpy)
S2cut = cutZero(S2numpy)
S3cut = cutZero(S3numpy)
plt.figure()
plt.plot(muCut,S1cut.imag,label=r'$S_1$',linewidth=2)
plt.plot(muCut,S2cut.imag,label=r'$S_2$',linewidth=2)
plt.plot(muCut,S3cut.imag,label=r'$S_3$',linewidth=2)
plt.xlim(-3,3)
plt.ylim(-6.5,6.5)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\mathrm{Im}(S)$')
plt.legend()
plt.tight_layout()
plt.show()