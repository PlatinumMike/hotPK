# -*- coding: utf-8 -*-
"""
Compute kernel on a structured 3D mesh for plotting. Only the part in parentheses, so excluding the prefactor in from.
You can actually do this far more efficiently, by exploiting the separation of variables.
So computing first on a radial, angular and axial grid, then taking the outer product.
But efficiency is not the goal here, it is just for making a plot, and possibly doing the FFT.
"""

import mpmath as mp
import numpy as np
import HfunctionGeneral as hfg
import Sfunctions as sf
import time
import matplotlib.pyplot as plt

def castToNumpy(mpArray, dim):
    if dim==1:
        return np.array(mpArray.tolist(), dtype='complex128')[0, :]
    else:
        return np.array(mpArray.tolist(), dtype='complex128')
    
def getPoints(resolutionX=100, resolutionY=100, resolutionZ=300, xMax=3, yMax=3, zMax=30, nMax=3, signCyc=1, vTOverC = 1.0e-3, freqRatio = 0.9, dataLayout = "normal", plotResults=False):
    """

    Parameters
    ----------
    resolutionX : int, optional
        mesh resolution in x direction. All of these resolutions must be even, because there is a singularity at the origin.
    resolutionY : TYPE, optional
        mesh resolution in y direction
    resolutionZ : TYPE, optional
        mesh resolution in z direction
    xMax : double, optional
        max of the range, already in normalized units: all distances divided by (2*thermal Larmor radius)
    yMax : double, optional
        max of y range
    zMax : double, optional
        max of z range
    nMax : int, optional
        max harmonic: -nMax, -nMax +1 ,...,0,1,2,..., nMax
    signCyc : int, optional
        sign of the cyclotron frequency (1 for positive particles, -1 for negative ions)
    vTOverC : double, optional
        thermal velocity divided by c, only needed for D3, D4, D5
    freqRatio : double, optional
        antenna frequency divided by cyclotron frequency
    dataLayout : string, optional
        if not using 'normal', it will wrap around the negative positions, which is useful for using the DFT.
    plotResults : bool, optional
        Make some basic plots of the kernel

    Returns
    -------
    x, y, z, kernel

    """
    #TODO: would be cleaner to return an object, then you can call plotting on it later on.

    cOverVT = 1.0/vTOverC
    NHarmonics = 2*nMax +  1

    #make mesh
    xRange = np.linspace(-xMax, xMax, resolutionX)
    yRange = np.linspace(-yMax, yMax, resolutionY)
    zRangeMP = mp.linspace(-zMax, zMax, resolutionZ)
    zRange = np.array(zRangeMP, dtype='float64')
    if dataLayout != "normal":
        xRange = np.fft.fftshift(xRange)
        yRange = np.fft.fftshift(yRange)
        zRange = np.fft.fftshift(zRange)
        zRangeMP = np.fft.fftshift(zRangeMP)
    
    #allocate space
    S1 = mp.matrix(NHarmonics, resolutionZ)
    S2 = mp.matrix(NHarmonics, resolutionZ)
    S3 = mp.matrix(NHarmonics, resolutionZ)
    HValues = np.zeros((7, NHarmonics, resolutionX, resolutionY), dtype='float64')
    radius = np.zeros((resolutionX, resolutionY), dtype='float64')
    angle = np.zeros((resolutionX, resolutionY), dtype='float64')
    kernel = np.zeros((4, 3, resolutionX, resolutionY, resolutionZ), dtype='complex128')

    startTime = time.time()
    print("Computing S functions")
    #compute kernel
    for iz in range(resolutionZ):
        #these are expensive to compute, but we can leverage the separation of variables, so they are the same for any x,y, thus do not need to be recomputed.
        for harmonic in range(-nMax, nMax+1):
            eps = np.sign(freqRatio-harmonic)*signCyc
            mu = zRangeMP[iz]*abs(freqRatio-harmonic)
            S1[harmonic+nMax, iz] = sf.getS1(eps, mu)
            S2[harmonic+nMax, iz] = sf.getS2(eps, mu)
            S3[harmonic+nMax, iz] = sf.getS3(eps, mu)
    
    S1numpy = castToNumpy(S1,2)
    S2numpy = castToNumpy(S2,2)
    S3numpy = castToNumpy(S3,2)
    print("Done, took ",time.time() - startTime," seconds")
    
    for ix in range(resolutionX):
        for iy in range(resolutionY):
            radius[ix, iy] = np.sqrt(xRange[ix]**2+yRange[iy]**2)
            angle[ix, iy] = np.arctan2(yRange[iy], xRange[ix])
    
    startTime = time.time()
    print("Computing H functions")
    #this can be done more efficiently, exploiting the fact that H only depends on radius, and you can use the mirror symmetry in terms of harmonic. But speed is not important for this simple example.
    for index in range(7):
        for harmonic in range(-nMax,nMax+1):
            for ix in range(resolutionX):
                for iy in range(resolutionY):
                    HValues[index, harmonic+nMax, ix, iy] = float(hfg.getH(radius[ix, iy], index, harmonic, False))
    print("Done, took ",time.time() - startTime," seconds")
    
    startTime = time.time()
    print("Filling up kernel")
    for ix in range(resolutionX):
        for iy in range(resolutionY):
            sine = np.sin(angle[ix, iy])
            cosine = np.cos(angle[ix, iy])
            sine2 = np.sin(2*angle[ix, iy])
            cosine2 = np.cos(2*angle[ix, iy])
            
            D3 = np.zeros((resolutionZ,), dtype='complex128')
            D4 = np.zeros((resolutionZ,), dtype='complex128')
            D5 = np.zeros((resolutionZ,), dtype='complex128')
            R0 = np.zeros((resolutionZ,), dtype='complex128')
            R0mod = np.zeros((resolutionZ,), dtype='complex128')
            R1 = np.zeros((resolutionZ,), dtype='complex128')
            R2 = np.zeros((resolutionZ,), dtype='complex128')
            R3 = np.zeros((resolutionZ,), dtype='complex128')
            R4 = np.zeros((resolutionZ,), dtype='complex128')
            R5 = np.zeros((resolutionZ,), dtype='complex128')
            for harmonicIndex in range(NHarmonics):
                D3 += -cOverVT*HValues[4,harmonicIndex,ix,iy] * S2numpy[harmonicIndex, :]
                D4 += 1.0j*np.sqrt(2)*signCyc*cOverVT*HValues[5,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                D5 += np.sqrt(2)*signCyc*cOverVT*HValues[6,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                R0 += 2*HValues[0,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                R0mod += -2*HValues[1,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                R1 += HValues[2,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                R2 += -1.0j*HValues[3,harmonicIndex,ix,iy] * S1numpy[harmonicIndex, :]
                R3 += -HValues[4,harmonicIndex,ix,iy] * S3numpy[harmonicIndex, :]
                R4 += -0.5j*np.sqrt(2)*signCyc*HValues[5,harmonicIndex,ix,iy] * S2numpy[harmonicIndex, :]
                R5 += -0.5*np.sqrt(2)*signCyc*HValues[6,harmonicIndex,ix,iy] * S2numpy[harmonicIndex, :]
    
            #row 1
            kernel[0, 0, ix, iy, :] = D4*cosine - D5*sine
            kernel[0, 1, ix, iy, :] = D4*sine + D5*cosine
            kernel[0, 2, ix, iy, :] = D3
            #row 2
            kernel[1, 0, ix, iy, :] = R1+0.5*R0-0.5*R0mod*cosine2
            kernel[1, 1, ix, iy, :] = R2 - 0.5*R0mod*sine2
            kernel[1, 2, ix, iy, :] = R4*cosine + R5*sine
            #row 3
            kernel[2, 0, ix, iy, :] = -R2 - 0.5*R0mod*sine2
            kernel[2, 1, ix, iy, :] = R1+0.5*R0+0.5*R0mod*cosine2
            kernel[2, 2, ix, iy, :] = R4*sine - R5*cosine
            #row 4
            kernel[3, 0, ix, iy, :] = R4*cosine - R5*sine
            kernel[3, 1, ix, iy, :] = R4*sine + R5*cosine
            kernel[3, 2, ix, iy, :] = R3
    
    print("Done, took ",time.time() - startTime," seconds")
    
    
    if plotResults:
        #some example plots, slicing the domain. Using the (4,3) element
        plt.figure()
        plt.plot(zRange,kernel[3,2,0,0,:].real)
        plt.xlabel(r"$\Delta z/(2\rho_T)$")
        
        plt.figure()
        plt.contourf(xRange,yRange,kernel[3,2,:,:,0].real.T,30)
        plt.xlabel(r"$\Delta x/(2\rho_T)$")
        plt.ylabel(r"$\Delta y/(2\rho_T)$")
        
        plt.figure()
        plt.contourf(xRange,zRange,kernel[3,2,:,0,:].real.T,30)
        plt.xlabel(r"$\Delta x/(2\rho_T)$")
        plt.ylabel(r"$\Delta z/(2\rho_T)$")
        
        plt.show()
    
    return xRange, yRange, zRange, kernel
