# -*- coding: utf-8 -*-
"""
Compute S functions.
Requires the mpmath package because it has the meijerg, but also because it
needs the arbitrary precision at large arguments.
"""

import numpy as np
import mpmath as mp

mp.mp.dps=100 #set precision

def getS1(eps,x):
    return mp.pi*1j*mp.meijerg([[],[]],[[0,0],[1/2]],x*x) + eps*mp.pi*(4*abs(x)*mp.hyper([],[3/2,3/2],x*x)-mp.sqrt(mp.pi)*mp.hyper([],[1/2,1],x*x))

def getS2(eps,x):
    return 2*mp.pi*eps*x*mp.meijerg([[],[]],[[0,0],[-1/2]],x*x) + mp.sign(x)*2j*mp.pi*(2*mp.sqrt(mp.pi)*abs(x)*mp.hyper([],[1,3/2],x*x)-mp.hyper([],[1/2,1/2],x*x))

def getS3(eps,x):
    return -2j*mp.pi*mp.meijerg([[],[]],[[0,1],[1/2]],x*x) + 4*mp.pi*eps*(abs(x)*mp.hyper([],[1/2,3/2],x*x)-mp.sqrt(mp.pi)*x*x*mp.hyper([],[3/2,2],x*x))
