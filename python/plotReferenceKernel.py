# -*- coding: utf-8 -*-
"""
Compute, and plot k-space kernel reference solution and FT of r-space kernel.
Taking the FT can be sped up considerably using the FFT.
However, in this script the FT is computed using Simpson's rule.
"""

import referenceKernel as rkern
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import generate3DPoints as g3dp
from scipy import integrate
import time

def manualDFT(function, x, y, z , kx, ky, kz):
    """
    Parameters
    ----------
    function : complex double
        values of f(x,y,z).
    x : doubles
        range of x.
    y : doubles
        range of y.
    z : doubles
        range of z.
    kx : double
        kx sample value.
    ky : double
        ky sample value.
    kz : double
        kz sample value.
    Returns
    -------
    complex
        F(kx,ky,kz) = FT(f)
    """
    
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    fLocal = function.copy()
    fLocal *= np.exp(-1.0j*kx*x[:,None,None])
    fLocal *= np.exp(-1.0j*ky*y[None,:,None])
    fLocal *= np.exp(-1.0j*kz*z[None,None,:])
    
    #Simpson rule
    fLocal = integrate.simpson(fLocal, dx=dz, axis=-1) #integrate over z
    fLocal = integrate.simpson(fLocal, dx=dy, axis=-1) #integrate over y
    F = integrate.simpson(fLocal, dx=dx) #integrate over x
    
    return F

def plotResultZ(row, col):
    plt.figure()
    plt.plot(kz,kernel[row,col,:].real, c='k', label='ref, Re')
    plt.plot(kz,kernel[row,col,:].imag, c='b', label='ref, Im')
    plt.plot(kz,kernel_FTMan[row,col,:].real, c='g', linestyle='--',label='num, Re')
    plt.plot(kz,kernel_FTMan[row,col,:].imag, c='r',linestyle='--',label='num, Im')
    plt.xlabel(r"$2\rho_T k_z$")
    plt.ylabel(r"$\tilde{{\sigma}}_{{{0}}}/(8\rho_T^3 \mathcal{{Q}})$".format(str(row+1)+str(col+1)))
    plt.legend()
    plt.tight_layout()
    
def plotResultX(row, col):
    plt.figure()
    plt.plot(kx2,kernel2[row,col,:].real, c='k', label='ref, Re')
    plt.plot(kx2,kernel2[row,col,:].imag, c='b', label='ref, Im')
    plt.plot(kx2,kernel_FTMan2[row,col,:].real, c='g', linestyle='--', label='num, Re')
    plt.plot(kx2,kernel_FTMan2[row,col,:].imag, c='r', linestyle='--', label='num, Im')
    plt.xlabel(r"$2\rho_T k_x$")
    plt.ylabel(r"$\tilde{{\sigma}}_{{{0}}}/(8\rho_T^3 \mathcal{{Q}})$".format(str(row+1)+str(col+1)))
    plt.legend()
    plt.tight_layout()

#physics inputs:
signCyc = 1
vTOverC = 0.001
freqRatio = 0.9

nMax = 3 #max harmonic

#mesh inputs:
resolutionX = 100  # make sure it is even to avoid singularity at the origin
# making it uniform is handy for if the FFT is used later on. Here I have not used FFT, but just simple numerical integration.
resolutionY = resolutionX
resolutionZ = 3000
xMax = 3
yMax = xMax
zMax = 300 #needs to be larger to resolve long wavelengths (small k), but that requires high precision in mpmath for stability of computing the S functions.

print("getting kernel on mesh, warning, this can take several hours")
xRange, yRange, zRange, kernelMesh = g3dp.getPoints(resolutionX, resolutionY, resolutionZ, xMax, yMax, zMax, nMax, signCyc, vTOverC, freqRatio, "normal", False)

#sample points k space
resolutionkX = 300
resolutionkZ = resolutionkX
kx = 0.1
ky = 0.15
kz = np.linspace(-4,4, resolutionkZ)

kx2 = np.linspace(-10,10, resolutionkX)
ky2 = -0.1
kz2 = 0.2

kernel = np.zeros((4,3,resolutionkZ),dtype='complex128')
kernel2 = np.zeros((4,3,resolutionkX),dtype='complex128')
kernel_FTMan = np.zeros_like(kernel)
kernel_FTMan2 = np.zeros_like(kernel2)

print("getting k-space kernel")
for sample in range(resolutionkZ):
    kernel[:,:,sample] = rkern.getKernel(kx, ky, kz[sample], nMax, signCyc, vTOverC, freqRatio)
for sample in range(resolutionkX):
    kernel2[:,:,sample] = rkern.getKernel(kx2[sample], ky2, kz2, nMax, signCyc, vTOverC, freqRatio)

print("Fourier transforming the kernel")
begin = time.time()
for row in range(4):
    for col in range(3):
        for ikz in range(resolutionkZ):
            kernel_FTMan[row,col,ikz] = manualDFT(kernelMesh[row,col,:,:,:], xRange, yRange, zRange, kx, ky, kz[ikz])
        for ikx in range(resolutionkX):
            kernel_FTMan2[row,col,ikx] = manualDFT(kernelMesh[row,col,:,:,:], xRange, yRange, zRange, kx2[ikx], ky2, kz2)
print("Done, took ",time.time()-begin," seconds")


plotResultZ(3,2)
plotResultX(3,2)

plt.show()