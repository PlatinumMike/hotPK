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
        reader = csv.reader(csvfile, delimiter=',')
        i=0
        for row in reader:
            data[i] = float(row[0])
            i+=1


def readCSVcomplex(fileName, data):
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        i=0
        for row in reader:
            data[i] = float(row[0])+1.0j*float(row[1])
            i+=1

runDir = "../build/"


#get number of lines
with open(runDir+"nodes.csv") as f:
    NR = sum(1 for line in f)


#load data
nodePositionR = np.zeros((NR,),'float64')
solutionVectorA = np.zeros((4*NR,),dtype='complex128')
readCSV(runDir+"nodes.csv", nodePositionR)
readCSVcomplex(runDir+"sol.csv", solutionVectorA)


#reshaping
speedOfLight = 299792458.0
omega = 40e6*2*np.pi #todo: get from inputs
pot = speedOfLight*solutionVectorA[:NR] #divided by c for better matrix conditioning, now retrieve potential so multiply times c.
AR = solutionVectorA[NR:2*NR]
Avarphi = solutionVectorA[2*NR:3*NR]
AZ = solutionVectorA[3*NR:]

EZ = 1.0j*omega*AZ
#todo: get ER,Evarphi, BR,BZ,Bvarphi...

#plotting
plt.figure()
plt.plot(nodePositionR,np.real(EZ),label='real part')
plt.plot(nodePositionR,np.imag(EZ),label='imaginary part')
plt.xlabel("R(m)")
plt.ylabel("EZ")
plt.legend()
plt.show()