#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results
"""

import numpy as np
import matplotlib.pyplot as plt
import csv


def readCSV(fileName, data):
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        i = 0
        for row in reader:
            data[i] = float(row[0])
            i += 1


def readCSVcomplex(fileName, data):
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        i = 0
        for row in reader:
            data[i] = float(row[0]) + 1.0j * float(row[1])
            i += 1


runDir = "../build/"


# get number of lines
with open(runDir + "nodes.csv") as f:
    NR = sum(1 for line in f)


# load data
nodePositionR = np.zeros((NR,), "float64")
solutionVectorA = np.zeros((4 * NR,), dtype="complex128")
readCSV(runDir + "nodes.csv", nodePositionR)
readCSVcomplex(runDir + "sol.csv", solutionVectorA)


# reshaping
speedOfLight = 299792458.0
omega = 51e6 * 2 * np.pi  # todo: get from inputs
# potential is divided by c for better matrix conditioning, now retrieve potential so multiply times c.
pot = speedOfLight * solutionVectorA[:NR]
AR = solutionVectorA[NR : 2 * NR]
Avarphi = solutionVectorA[2 * NR : 3 * NR]
AZ = solutionVectorA[3 * NR :]

EZ = 1.0j * omega * AZ
# todo: get ER,Evarphi, BR,BZ,Bvarphi...

# plotting
plt.figure()
plt.plot(nodePositionR, np.real(pot) / speedOfLight, label="real part")
plt.plot(nodePositionR, np.imag(pot) / speedOfLight, label="imaginary part")
plt.xlabel("R(m)")
plt.ylabel(r"$\Phi$ (Vs/m)")
plt.legend()

plt.figure()
plt.plot(nodePositionR, np.real(AZ), label="real part")
plt.plot(nodePositionR, np.imag(AZ), label="imaginary part")
plt.xlabel("R(m)")
plt.ylabel(r"$A_Z$ (Vs/m)")
plt.legend()
plt.show()
