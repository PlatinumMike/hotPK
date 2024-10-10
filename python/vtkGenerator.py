# -*- coding: utf-8 -*-
"""
Generate and save data for Paraview.
Using the pyevtk (export to VTK) package to save some datasets in a more Paraview friendly format.
"""

import generate3DPoints as g3dp
import numpy as np
from pyevtk.hl import gridToVTK


def saveData(fileName: str, x, y, z, kernel):
    print("saving vtk data now")
    # data needs to be a double, and contiguous, so storing it in two new arrays.
    kernReal = np.empty(kernel.shape, dtype="float64")
    kernReal = kernel.real.copy()
    kernImag = np.empty(kernel.shape, dtype="float64")
    kernImag = kernel.imag.copy()
    gridToVTK(
        fileName,
        x,
        y,
        z,
        pointData={
            "kernel11_re": kernReal[0, 0, :, :, :],
            "kernel11_im": kernImag[0, 0, :, :, :],
            "kernel12_re": kernReal[0, 1, :, :, :],
            "kernel12_im": kernImag[0, 1, :, :, :],
            "kernel13_re": kernReal[0, 2, :, :, :],
            "kernel13_im": kernImag[0, 2, :, :, :],
            "kernel21_re": kernReal[1, 0, :, :, :],
            "kernel21_im": kernImag[1, 0, :, :, :],
            "kernel22_re": kernReal[1, 1, :, :, :],
            "kernel22_im": kernImag[1, 1, :, :, :],
            "kernel23_re": kernReal[1, 2, :, :, :],
            "kernel23_im": kernImag[1, 2, :, :, :],
            "kernel31_re": kernReal[2, 0, :, :, :],
            "kernel31_im": kernImag[2, 0, :, :, :],
            "kernel32_re": kernReal[2, 1, :, :, :],
            "kernel32_im": kernImag[2, 1, :, :, :],
            "kernel33_re": kernReal[2, 2, :, :, :],
            "kernel33_im": kernImag[2, 2, :, :, :],
            "kernel41_re": kernReal[3, 0, :, :, :],
            "kernel41_im": kernImag[3, 0, :, :, :],
            "kernel42_re": kernReal[3, 1, :, :, :],
            "kernel42_im": kernImag[3, 1, :, :, :],
            "kernel43_re": kernReal[3, 2, :, :, :],
            "kernel43_im": kernImag[3, 2, :, :, :],
        },
    )
    print("done saving")


# mesh inputs:
resolutionX = 100  # make sure it is even to avoid singularity at the origin
# making it uniform is handy for taking the FFT later on.
resolutionY = resolutionX
resolutionZ = 300
xMax = 3
yMax = xMax
zMax = 30  # longer dispersion along the z direction, make larger than xMax, yMax
nMax = 3  # try larger values as well to check for convergence.

# physics inputs:
signCyc = 1  # I will use ions here, so positive
vTOverC = 0.001  # v_T / c
freqRatio = 0.9  # \omega/\Omega.

xRange, yRange, zRange, kernel = g3dp.getPoints(
    resolutionX, resolutionY, resolutionZ, xMax, yMax, zMax, nMax, signCyc, vTOverC, freqRatio, "normal", True
)

# saving data now for further analysis in Paraview
saveData("kernelData2", xRange, yRange, zRange, kernel)
