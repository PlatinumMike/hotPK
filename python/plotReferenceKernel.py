# -*- coding: utf-8 -*-
"""
Compute configuration space kernel on a grid, then FT and compare against known results.
Since it is a sum over harmonics, and the kernel is separable, this can be optimised.
So the transform over z direction is independent of that over x,y.
Note that this is just a validation, speed is not the goal.

You have to use a non-uniform grid in z to accurately resolve the spike near the 
origin with a reasonable number of grid points.
This rules out the use of the FFT (needs uniform grid).
But perhaps you could use some coordinate transformation in which the singularity
in S1 is removed. E.g. switching to spherical coordinates tacks on a factor r^2 (jacobian).
This kills off the singularities in both the S and H functions!
Here polar coordinates in x,y, are used to avoid the singularity.
Also exploiting the fact that the separation of values gives you the ability to 
pull out the S functions, FT them, and FT the left over bit. Then merge together.

Btw, if the used resolution or max of the range is insufficient you will get 
some typical problems that fall into two classes: (1) the function is the same as
the reference value, but shifted by some constant amount, in this case the peak 
is not well resolved, use more mesh points (near the peak). (2) the shape does
not match the reference value, in this case the range max is insufficient, and
thus some important non-zero part of the S functions is cut out. Increase zMax.
"""

import numpy as np
import referenceKernel as rkern
import matplotlib.pyplot as plt
import HfunctionGeneral as hfg
import Sfunctions as sf
import mpmath as mp
from scipy import integrate
import time

plt.rcParams.update({"font.size": 14})


def castToNumpy(mpArray, dim: int):
    if dim == 1:
        return np.array(mpArray.tolist(), dtype="complex128")[0, :]
    else:
        return np.array(mpArray.tolist(), dtype="complex128")


def getS(z: np.ndarray, harmonic: int, freqRatio: float, version: int):
    S = mp.matrix(1, z.size)
    for iz in range(z.size):
        eps = np.sign(freqRatio - harmonic) * signCyc
        mu = z[iz] * abs(freqRatio - harmonic)
        if version == 1:
            S[0, iz] = sf.getS1(eps, mu)
        elif version == 2:
            S[0, iz] = sf.getS2(eps, mu)
        else:
            S[0, iz] = sf.getS3(eps, mu)
    Snumpy = castToNumpy(S, 1)
    return Snumpy


def getKernel(HValues: np.ndarray, angle: float):
    sine = np.sin(angle)
    cosine = np.cos(angle)
    sine2 = np.sin(2 * angle)
    cosine2 = np.cos(2 * angle)

    D3 = -cOverVT * HValues[4]
    D4 = 1.0j * np.sqrt(2) * signCyc * cOverVT * HValues[5]
    D5 = np.sqrt(2) * signCyc * cOverVT * HValues[6]
    R0 = 2 * HValues[0]
    R0mod = -2 * HValues[1]
    R1 = HValues[2]
    R2 = -1.0j * HValues[3]
    R3 = -HValues[4]
    R4 = -0.5j * np.sqrt(2) * signCyc * HValues[5]
    R5 = -0.5 * np.sqrt(2) * signCyc * HValues[6]

    kernel = np.zeros((4, 3), dtype="complex128")

    # row 1
    kernel[0, 0] = D4 * cosine - D5 * sine
    kernel[0, 1] = D4 * sine + D5 * cosine
    kernel[0, 2] = D3
    # row 2
    kernel[1, 0] = R1 + 0.5 * R0 - 0.5 * R0mod * cosine2
    kernel[1, 1] = R2 - 0.5 * R0mod * sine2
    kernel[1, 2] = R4 * cosine + R5 * sine
    # row 3
    kernel[2, 0] = -R2 - 0.5 * R0mod * sine2
    kernel[2, 1] = R1 + 0.5 * R0 + 0.5 * R0mod * cosine2
    kernel[2, 2] = R4 * sine - R5 * cosine
    # row 4
    kernel[3, 0] = R4 * cosine - R5 * sine
    kernel[3, 1] = R4 * sine + R5 * cosine
    kernel[3, 2] = R3
    return kernel


def fillZValues(S1, S2, S3):
    kernel = np.zeros((4, 3), dtype="complex128")

    kernel[0, 0] = S1  # D4, D5 share S1
    kernel[0, 1] = S1
    kernel[0, 2] = S2  # D3 has S2

    kernel[1, 0] = S1  # R0, R0mod, R1 have S1
    kernel[1, 1] = S1  # etc...
    kernel[1, 2] = S2

    kernel[2, 0] = S1
    kernel[2, 1] = S1
    kernel[2, 2] = S2

    kernel[3, 0] = S2
    kernel[3, 1] = S2
    kernel[3, 2] = S3  # R3

    return kernel


def manualFT(function, z, kz):
    """
    Parameters
    ----------
    function : complex double
        values of f(z)
    z : doubles
        range of z.
    kz : double
        kz sample value.
    Returns
    -------
    complex
        F(kz) = FT(f)
    """
    return integrate.simpson(y=function * np.exp(-1.0j * z * kz), x=z)


def manualFT2(f, rRange, angles, kx: float, ky: float):
    """

    Parameters
    ----------
    f : complex doubles
        values of f(xi, alpha)
    rRange : doubles
        range of xi values.
    angles : doubles
        range of angles.
    kx : double
        kx sample value.
    ky : double
        ky sample value.

    Returns
    -------
    complex
        F(kx,ky) = FT(f)
    """

    kperp = np.sqrt(kx**2 + ky**2)
    psi = np.arctan2(ky, kx)
    integrand = f * np.exp(-1.0j * kperp * rRange[:, None] * np.cos(angles[None, :] - psi))
    # integrate over angle, then integrate over radial position.
    return integrate.simpson(y=integrate.simpson(y=integrand, x=angles, axis=-1), x=rRange)


def plotResultZ(row: int, col: int):
    plt.figure()
    plt.plot(kz, kernel[row, col, :].real, c="k", label="ref, Re")
    plt.plot(kz, kernel[row, col, :].imag, c="b", label="ref, Im")
    plt.plot(kz, kernel_FT[row, col, :].real, c="g", linestyle="--", label="num, Re")
    plt.plot(kz, kernel_FT[row, col, :].imag, c="r", linestyle="--", label="num, Im")
    plt.xlabel(r"$2\rho_T k_z$")
    plt.ylabel(r"$\tilde{{\sigma}}_{{{0}}}/(8\rho_T^3 \mathcal{{Q}})$".format(str(row + 1) + str(col + 1)))
    plt.legend()
    plt.tight_layout()


def plotResultX(row: int, col: int):
    plt.figure()
    plt.plot(kx2, kernel2[row, col, :].real, c="k", label="ref, Re")
    plt.plot(kx2, kernel2[row, col, :].imag, c="b", label="ref, Im")
    plt.plot(kx2, kernel_FT2[row, col, :].real, c="g", linestyle="--", label="num, Re")
    plt.plot(kx2, kernel_FT2[row, col, :].imag, c="r", linestyle="--", label="num, Im")
    plt.xlabel(r"$2\rho_T k_x$")
    plt.ylabel(r"$\tilde{{\sigma}}_{{{0}}}/(8\rho_T^3 \mathcal{{Q}})$".format(str(row + 1) + str(col + 1)))
    plt.legend()
    plt.tight_layout()


# physics inputs:
vTOverC = 0.001
freqRatio = 0.9
signCyc = np.sign(freqRatio)  # defining antenna freq to be postive.
cOverVT = 1.0 / vTOverC

nMin = -3  # min harmonic
nMax = 3  # max harmonic
NHarmonics = nMax - nMin + 1

# k mesh
resKz = 1000
kzMax = 4
# already normalized, so times 2 \rho_T, and distances z are divided by 2\rho_T.
kz = np.linspace(-kzMax, kzMax, resKz)
minWavelengthZ = 2 * np.pi / kzMax
kx = 0.1
ky = 0.15

# second example, but now varying kx
resKx = resKz
kxMax = 10
kx2 = np.linspace(-kxMax, kxMax, resKx)
ky2 = -0.1
kz2 = 0.2
minWavelengthX = 2 * np.pi / kxMax

# mesh inputs
numPointsPerWavelength = 30
rMax = 3
numAngles = 100
resR = int(numPointsPerWavelength * 2 * rMax / minWavelengthX)
print("Using rMax = {}, resR = {}".format(rMax, resR))
angles = np.linspace(0, 2 * np.pi, numAngles)
rRange = np.linspace(0, rMax, resR)  # 0 is not a problem now, because of Jacobian.

kernel = np.zeros((4, 3, resKz), dtype="complex128")  # reference values
kernel2 = np.zeros((4, 3, resKx), dtype="complex128")

kernel_FT = np.zeros_like(kernel)  # test values
kernel_FT2 = np.zeros_like(kernel2)

S1_FT = np.zeros((resKz,), dtype="complex128")  # FT of the S functions
S2_FT = np.zeros_like(S1_FT)
S3_FT = np.zeros_like(S1_FT)
S1_FT2 = 0j  # second experiment just has one kz value
S2_FT2 = 0j
S3_FT2 = 0j
zValues = np.zeros((4, 3, resKz), dtype="complex128")
zValues2 = np.zeros((4, 3), dtype="complex128")

# XY grid values
rValues = np.zeros((4, 3, resR, numAngles), dtype="complex128")
HValues = np.zeros((7, resR), dtype="float64")

print("getting k-space kernel")
for sample in range(resKz):
    kernel[:, :, sample] = rkern.getKernel(kx, ky, kz[sample], nMin, nMax, signCyc, vTOverC, freqRatio)
for sample in range(resKx):
    kernel2[:, :, sample] = rkern.getKernel(kx2[sample], ky2, kz2, nMin, nMax, signCyc, vTOverC, freqRatio)

# z mesh, this 'z' is already divided by 2\rho_T
numPointsPerWavelength = 100
muMax = 12
# resolution is based on the \mu dependency, not z. But the resolution to resolve a specific kz depends on the z resolution. So needs to be set for each harmonic separately, in order to get an accurate result.

print("getting kernel on z mesh now, and Fourier transforming, warning, this can take several minutes")
for harmonicIndex in range(NHarmonics):
    now = time.time()
    harmonic = nMin + harmonicIndex
    zMax = muMax / np.abs(freqRatio - harmonic)
    resZ = int(numPointsPerWavelength * 2 * zMax / minWavelengthZ)
    resZ += resZ % 2  # must be even, so that it can be split into two pieces.
    print("Harmonic = {}, using zMax = {}, resZ = {}".format(harmonic, zMax, resZ))
    peaking = 3.0
    tmp = np.linspace(0, 1, resZ // 2 + 1)[1:] ** peaking
    # has no zero, but peaked near the origin to resolve spike better.
    zRange = np.concatenate((-np.flip(tmp), tmp)) * zMax
    Sfunc1 = getS(zRange, harmonic, freqRatio, 1)  # compute S functions
    Sfunc2 = getS(zRange, harmonic, freqRatio, 2)
    Sfunc3 = getS(zRange, harmonic, freqRatio, 3)

    # transform in z direction
    for j in range(resKz):
        S1_FT[j] = manualFT(Sfunc1, zRange, kz[j])
        S2_FT[j] = manualFT(Sfunc2, zRange, kz[j])
        S3_FT[j] = manualFT(Sfunc3, zRange, kz[j])
        zValues[:, :, j] = fillZValues(S1_FT[j], S2_FT[j], S3_FT[j])
    S1_FT2 = manualFT(Sfunc1, zRange, kz2)
    S2_FT2 = manualFT(Sfunc2, zRange, kz2)
    S3_FT2 = manualFT(Sfunc3, zRange, kz2)
    zValues2 = fillZValues(S1_FT2, S2_FT2, S3_FT2)

    for i in range(7):
        for j in range(resR):
            # we might as well already include factor xi already here.
            HValues[i, j] = hfg.getH(rRange[j], i, harmonic, jacobian=True)
    print("Getting values in the XY plane now")
    for i in range(resR):
        for j in range(numAngles):
            rValues[:, :, i, j] = getKernel(HValues[:, i], angles[j])

    # transform in x,y plane
    kernelXY_FT = np.zeros((4, 3), dtype="complex128")
    kernelXY_FT2 = np.zeros((4, 3, resKx), dtype="complex128")
    for i in range(4):
        for j in range(3):
            kernelXY_FT[i, j] = manualFT2(rValues[i, j, :, :], rRange, angles, kx, ky)
            for k in range(resKx):
                kernelXY_FT2[i, j, k] = manualFT2(rValues[i, j, :, :], rRange, angles, kx2[k], ky2)

    # reconstruct kernel, sum over every harmonic
    for i in range(4):
        for j in range(3):
            for k in range(resKz):
                kernel_FT[i, j, k] += kernelXY_FT[i, j] * zValues[i, j, k]
            for l in range(resKx):
                kernel_FT2[i, j, l] += kernelXY_FT2[i, j, l] * zValues2[i, j]
    print("Done, took {:.1f} s".format(time.time() - now))


# plotting
plotResultZ(3, 2)
plotResultX(3, 2)

plt.show()
