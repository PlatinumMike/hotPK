# -*- coding: utf-8 -*-
"""
Compute reference kernel values, in k-space.
Divided by (8*\rho_T^3*\mathcal{Q})
"""

import plasmapy.dispersion as dp
from scipy import special as sp
import numpy as np

def getKernel(kx, ky, kz, nMin=-3, nMax=3, signCyc=1, vToverC=1.0e-3, freqRatio=0.9):
    """

    Parameters
    ----------
    kx : double
        wave number kx, multiplied by 2 thermal larmor radii.
    ky : double
        wave number ky, multiplied by 2 thermal larmor radii.
    kz : double
        wave number kz, multiplied by 2 thermal larmor radii.
    nMax: int
        max harmonic
    physics inputs below. There is more, but these are the only FREE parameters to compute whatever is inside of the parentheses, if you normalize to kill off the prefactor.
    Actually, you may assume the antenna freq is positive, so then the sign of cyc can be derived from freqRatio.
    signCyc : int, optional
        sign of cyclotron frequency.
    vToverC : double, optional
        thermal velocity divided by speed of light
    freqRatio : double, optional
        antenna frequency divided by cyclotron freq.
        
    Warning: the edge cases where kz=0, or kx=ky=0 are not handled

    Returns
    -------
    kernel : complex
   
    """
    cOverVT = 1/vToverC
    kperp = np.sqrt(kx**2 + ky**2)
    FLR_parameter = 0.125*kperp**2
    exp = np.exp(-FLR_parameter)
    kzAbs = abs(kz)
    signKz = np.sign(kz)
    psi = np.arctan2(ky,kx)
    cos = np.cos(psi)
    sin = np.sin(psi)
    
    kernel = np.zeros((4,3),dtype='complex128')
    
    C3 = 0.0j
    C4 = 0.0j
    C5 = 0.0j
    K0 = 0.0j
    K1 = 0.0j
    K2 = 0.0j
    K3 = 0.0j
    K4 = 0.0j
    K5 = 0.0j
    for n in range(nMin,nMax+1):
        zeta = 2*signCyc*(freqRatio-n)/kzAbs
        PDF = dp.plasma_dispersion_func(zeta)
        dPDF = dp.plasma_dispersion_func_deriv(zeta)
        #this will fail if kperp becomes really large, e.g. 10000, cause then exp underflows and Iv overflows, you could use exponentially scaled bessel functions.
        besselI = sp.iv(n, FLR_parameter)
        besselIdiff = besselI - sp.ivp(n, FLR_parameter)
        
        C3 += -cOverVT*np.pi*exp/(4*kzAbs)*besselI*signKz*dPDF
        C4 += signCyc*cOverVT*np.pi*kperp*exp / (8*kzAbs)*n*besselI/FLR_parameter*PDF
        C5 += -1.0j*signCyc*cOverVT*np.pi*kperp*exp / (8*kzAbs)*besselIdiff*PDF
        K0 += np.pi*exp/(2*kzAbs)*FLR_parameter*besselIdiff*PDF
        K1 += np.pi*exp/(4*kzAbs)*n**2*besselI/FLR_parameter*PDF
        K2 += -1.0j*np.pi*exp/(4*kzAbs)*n*besselIdiff*PDF
        K3 += -np.pi*exp/(4*kzAbs)*besselI*zeta*dPDF
        K4 += -signCyc*np.pi*kperp*exp/(16*kzAbs)*n*besselI/FLR_parameter*signKz*dPDF
        K5 += 1.0j*signCyc*np.pi*kperp*exp/(16*kzAbs)*besselIdiff*signKz*dPDF
    
    #row 1
    kernel[0, 0] = C4*cos - C5*sin
    kernel[0, 1] = C4*sin + C5*cos
    kernel[0, 2] = C3
    #row 2
    kernel[1, 0] = K1 + K0*sin**2
    kernel[1, 1] = K2-K0*cos*sin
    kernel[1, 2] = K4*cos+K5*sin
    #row 3
    kernel[2, 0] = -K2 - K0*cos*sin
    kernel[2, 1] = K1 + K0*cos**2
    kernel[2, 2] = K4*sin - K5*cos
    #row 4
    kernel[3, 0] = K4*cos - K5*sin
    kernel[3, 1] = K4*sin + K5*cos
    kernel[3, 2] = K3

    return kernel
