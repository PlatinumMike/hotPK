# -*- coding: utf-8 -*-
"""
Contains H functions and their integrals
"""

import numpy as np
from scipy import special

a0H = [-1,-1,1 ,1 ,-1,0 ,0]
a1H = [0 , 0, 0,0 ,0 ,1 ,1]
b0H = np.array([[1/8,1/8,1/8,1/8],
               [1/8,1/8,1/8,1/8],
               [0,-1,-6,-15],
               [0,-1,-2,-3],
               [1/4,1/4,1/4,1/4],
               [0,1/np.sqrt(2),np.sqrt(2),3/np.sqrt(2)],
               [1/(2*np.sqrt(2)),1/(2*np.sqrt(2)),1/(2*np.sqrt(2)),1/(2*np.sqrt(2))]])/np.sqrt(np.pi)
b1H = np.array([[-1/4,1/4,15/4,41/4],
               [1/4,-1/4,-7/4,-17/4],
               [0,0,-4,-32],
               [0,0,-4,-20],
               [0,0,2,6],
               [0,0,2*np.sqrt(2),10*np.sqrt(2)],
               [0,0,2*np.sqrt(2),6*np.sqrt(2)]])/np.sqrt(np.pi)
b2H = np.array([[0,0,0,12],
               [0,0,0,-8],
               [0,0,0,-8],
               [0,0,0,-8],
               [0,0,0,4],
               [0,0,0,4*np.sqrt(2)],
               [0,0,0,4*np.sqrt(2)]])/np.sqrt(np.pi)
c0H = np.array([[0,-1/2,-1,-3/2],
               [0,0,0,0],
               [0,1/2,1,3/2],
               [0,0,0,0],
               [0,-1/2,-1,-3/2],
               [0,-1/np.sqrt(2),-2*np.sqrt(2),-9/np.sqrt(2)],
               [0,-1/np.sqrt(2),-np.sqrt(2),-3/np.sqrt(2)]])
c1H = np.array([[0,0,-4,-16],[0,0,2,8],[0,1,8,27],[0,1,4,9],
               [0,0,-2,-8],[0,0,-2*np.sqrt(2),-12*np.sqrt(2)],[0,0,-2*np.sqrt(2),-8*np.sqrt(2)]])
c2H = np.array([[0,0,0,-12],[0,0,0,8],[0,0,4,36],[0,0,4,24],[0,0,0,-4],
               [0,0,0,-4*np.sqrt(2)],[0,0,0,-4*np.sqrt(2)]])
c3H = np.array([[0,0,0,0],
               [0,0,0,0],
               [0,0,0,8],
               [0,0,0,8],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0]])

a0U = [1, 1, 1 ,1 ,1,0 , 0]
a1U = [2 , 0, 0, 0 ,0 ,3 ,3]
b0U = np.array([[1/8,1/8,1/8,1/8],
               [-1/8,1/8,1/8,1/8],
               [0,-1/8,0,0],
               [0,1/8,0,0],
               [0,1/4,1/4,1/4],
               [0,-1/(6*np.sqrt(2)),-1/(15*np.sqrt(2)),-3/(70*np.sqrt(2))],
               [-1/(4*np.sqrt(2)),1/(12*np.sqrt(2)),1/(60*np.sqrt(2)),1/(140*np.sqrt(2))]])/np.sqrt(np.pi)
b1U = np.array([[0,0,1,3],
               [0,0,-1/2,-4/3],
               [0,-1/4,-5/3,-9/2],
               [0,-1/4,-2/3,-1],
               [0,0,1/2,5/3],
               [0,1/(3*np.sqrt(2)),7*np.sqrt(2)/15,51/(35*np.sqrt(2))],
               [0,1/(3*np.sqrt(2)),2*np.sqrt(2)/15,9/(35*np.sqrt(2))]])/np.sqrt(np.pi)
b2U = np.array([[0,0,0,2],
               [0,0,0,-4/3],
               [0,0,-2/3,-11/2],
               [0,0,-2/3,-7/2],
               [0,0,0,2/3],
               [0,0,2*np.sqrt(2)/5,74*np.sqrt(2)/35],
               [0,0,2*np.sqrt(2)/5,46*np.sqrt(2)/35]])/np.sqrt(np.pi)
b3U = np.array([[0,0,0,0],
               [0,0,0,0],
               [0,0,0,-1],
               [0,0,0,-1],
               [0,0,0,0],
               [0,0,0,4*np.sqrt(2)/7],
               [0,0,0,4*np.sqrt(2)/7]])/np.sqrt(np.pi)

c0U = np.array([[0,-1/4,-1/2,-3/4],
               [-1/8,0,0,0],
               [0,-1/16,0,0],
               [0,1/16,0,0],
               [-1/8,0,0,0],
               [0,-1/(3*np.sqrt(2)),-2*np.sqrt(2)/3,-3/np.sqrt(2)],
               [0,-1/(3*np.sqrt(2)),-np.sqrt(2)/3,-1/np.sqrt(2)]])
c1U = np.array([[0,0,-1,-4],
               [0,0,0,0],
               [0,1/4,1/2,3/4],
               [0,0,0,0],
               [0,-1/4,-1/2,-3/4],
               [0,0,-2*np.sqrt(2)/5,-12*np.sqrt(2)/5],
               [0,0,-2*np.sqrt(2)/5,-8*np.sqrt(2)/5]])
c2U = np.array([[0,0,0,-2],
               [0,0,1/2,2],
               [0,1/4,2,27/4],
               [0,1/4,1,9/4],
               [0,0,-1/2,-2],
               [0,0,0,-4*np.sqrt(2)/7],
               [0,0,0,-4*np.sqrt(2)/7]])
c3U = np.array([[0,0,0,0],
               [0,0,0,4/3],
               [0,0,2/3,6],
               [0,0,2/3,4],
               [0,0,0,-2/3],
               [0,0,0,0],
               [0,0,0,0]])
c4U = np.array([[0,0,0,0],
               [0,0,0,0],
               [0,0,0,1],
               [0,0,0,1],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0]])

def getH(x,index,n,jacobian=False):        
    '''
    Parameters
    ----------
    x : position
    index : select H0, H0breve, H1, H2, H3, H4, H5
    n : harmonic, only mapped in the range -3 to 3
    jacobian : bool, optional
        if True, multiply by x so that the singularity is avoided.
        The default is False.

    Returns
    -------
    double
    '''
    
    
    signFlip = 1
    if n < 0 and (index == 3 or index == 5):
        signFlip = -1
        
    n = abs(n) #H_-n = H_n, except for H_2, H_4!
    
    a0local = a0H[index]
    a1local = a1H[index]
    if jacobian:
        a0local+=1
        a1local+=1
    
    poly1 = b0H[index,n] + b1H[index,n]*x*x + b2H[index,n]*x**4
    poly2 = c0H[index,n] + c1H[index,n]*x*x + c2H[index,n]*x**4 + c3H[index,n]*x**6
    
    return signFlip * (x**a0local * poly1 * np.exp(-x*x) + x**a1local * poly2 * special.erfc(x))

def getU(x,index,n):        
    '''
    U = \int dx x H(x)
    
    Parameters
    ----------
    x : position
    index : select U0, U0breve, U1, U2, U3, U4, U5
    n : harmonic, only mapped in the range -3 to 3

    Returns
    -------
    double
    '''
    
    signFlip = 1
    if n < 0 and (index == 3 or index == 5):
        signFlip = -1
        
    n = abs(n) #H_-n = H_n, except for H_2, H_4!
    
    a0local = a0U[index]
    a1local = a1U[index]
    
    poly1 = b0U[index,n] + b1U[index,n]*x*x + b2U[index,n]*x**4 + b3U[index,n]*x**6
    poly2 = c0U[index,n] + c1U[index,n]*x*x + c2U[index,n]*x**4 + c3U[index,n]*x**6 + c4U[index,n]*x**8
    
    return signFlip * (x**a0local * poly1 * np.exp(-x*x) + x**a1local * poly2 * special.erfc(x))

def getV(x,index,n):        
    '''
    V = \int dx x^2 H(x)
    
    Parameters
    ----------
    x : position
    index : select V0, V0breve, V1, V2, V3, V4, V5
    n : harmonic, only mapped in the range -3 to 3

    Returns
    -------
    double
    '''
    
    
    signFlip = 1
    if n < 0 and (index == 3 or index == 5):
        signFlip = -1
        
    n = abs(n) #H_-n = H_n, except for H_2, H_4!
    
    #copy-paste from Mathematica's CForm[]. Not per se optimized for rapid evaluation, but that can be improved later on.
    if index == 0:
        if n==0:
            return (1 + 2*x*x)/(16.*np.sqrt(np.pi))*np.exp(-x*x)
        elif n==1:
            return (-1 + 2*x*x)/(48.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*special.erfc(x))/6.
        elif n==2:
            return (-1 + 14*x*x + 192*pow(x,4))/(240.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*(5 + 12*x*x)*special.erfc(x))/15.
        elif n==3:
            return (-1 + 34*x*x + 1312*pow(x,4) + 960*pow(x,6))/(560.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*(35 + 224*x*x + 120*pow(x,4))*special.erfc(x))/70.
    elif index == 1:
        if n==0:
            return -0.0625*(3 + 2*x*x)/(np.sqrt(np.pi))*np.exp(-x*x)
        elif n==1:
            return (1 + 2*x*x)/(16.*np.sqrt(np.pi))*np.exp(-x*x)
        elif n==2:
            return -0.0125*(-1 - 6*x*x + 32*pow(x,4))/(np.sqrt(np.pi))*np.exp(-x*x) + (2*pow(x,5)*special.erfc(x))/5.
        elif n==3:
            return -0.0017857142857142857*(-3 - 38*x*x + 576*pow(x,4) + 640*pow(x,6))/(np.sqrt(np.pi))*np.exp(-x*x) + (8*pow(x,5)*(7 + 5*x*x)*special.erfc(x))/35.
    elif index == 2:
        if n==0:
            return 0
        elif n==1:
            return -0.06666666666666667*(1 + x*x + 3*pow(x,4))/(np.sqrt(np.pi))*np.exp(-x*x) + (pow(x,3)*(5 + 6*x*x)*special.erfc(x))/30.
        elif n==2:
            return (-2*(-2 - 2*x*x + 69*pow(x,4) + 30*pow(x,6)))/(105.*np.sqrt(np.pi))*np.exp(-x*x) + (pow(x,3)*(35 + 168*x*x + 60*pow(x,4))*special.erfc(x))/105.
        elif n==3:
            return -0.0031746031746031746*(-3 - 3*x*x + 1101*pow(x,4) + 1480*pow(x,6) + 280*pow(x,8))/(np.sqrt(np.pi))*np.exp(-x*x) + (pow(x,3)*(315 + 3402*x*x + 3240*pow(x,4) + 560*pow(x,6))*special.erfc(x))/630.
    elif index == 3:
        if n==0:
            return 0
        elif n==1:
            return signFlip*(-0.1*((-1 + x)*(1 + x)*(1 + 2*x*x))/(np.sqrt(np.pi))*np.exp(-x*x) + (pow(x,5)*special.erfc(x))/5.)
        elif n==2:
            return signFlip*(-0.02857142857142857*(1 + x*x + 18*pow(x,4) + 20*pow(x,6))/(np.sqrt(np.pi))*np.exp(-x*x) + (4*pow(x,5)*(7 + 5*x*x)*special.erfc(x))/35.)
        elif n==3:
            return signFlip*(-0.0015873015873015873*(3 + 3*x*x + 474*pow(x,4) + 1880*pow(x,6) + 560*pow(x,8))/(np.sqrt(np.pi))*np.exp(-x*x) + (pow(x,5)*(567 + 1080*x*x + 280*pow(x,4))*special.erfc(x))/315.)
    elif index == 4:
        if n==0:
            return -0.125*1/(np.sqrt(np.pi))*np.exp(-x*x)
        elif n==1:
            return (1 + 4*x*x)/(24.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*special.erfc(x))/6.
        elif n==2:
            return ((1 + 4*x*x)*(1 + 12*x*x))/(120.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*(5 + 6*x*x)*special.erfc(x))/15.
        elif n==3:
            return (1 + 36*x*x + 368*pow(x,4) + 160*pow(x,6))/(280.*np.sqrt(np.pi))*np.exp(-x*x) - (pow(x,3)*(35 + 112*x*x + 40*pow(x,4))*special.erfc(x))/70.
    elif index == 5:
        if n==0:
            return 0
        elif n==1:
            return signFlip*((x*(-1 + 2*x*x))/(8.*np.sqrt(2*np.pi))*np.exp(-x*x) - ((1 + 4*pow(x,4))*special.erfc(x))/(16.*np.sqrt(2)))
        elif n==2:
            return signFlip*((np.sqrt(2/np.pi)*pow(x,3)*(1 + x*x))/3*np.exp(-x*x) - (pow(x,4)*(3 + 2*x*x)*special.erfc(x))/(3.*np.sqrt(2)))
        elif n==3:
            return signFlip*((pow(x,3)*(2 + 7*x*x + 2*pow(x,4)))/(2.*np.sqrt(2*np.pi))*np.exp(-x*x) - (pow(x,4)*(9 + 16*x*x + 4*pow(x,4))*special.erfc(x))/(4.*np.sqrt(2)))
    elif index == 6:
        if n==0:
            return -0.25*x/np.sqrt(2*np.pi)*np.exp(-x*x) - special.erfc(x)/(8.*np.sqrt(2))
        elif n==1:
            return (x*(1 + 2*x*x))/(8.*np.sqrt(2*np.pi))*np.exp(-x*x) - ((-1 + 2*x*x)*(1 + 2*x*x)*special.erfc(x))/(16.*np.sqrt(2))
        elif n==2:
            return (pow(x,3)*(1 + 4*x*x))/(6.*np.sqrt(2*np.pi))*np.exp(-x*x) - (pow(x,4)*(3 + 4*x*x)*special.erfc(x))/(6.*np.sqrt(2))
        elif n==3:
            return (pow(x,3)*(1 + 13*x*x + 6*pow(x,4)))/(6.*np.sqrt(2*np.pi))*np.exp(-x*x) - (pow(x,4)*(9 + 32*x*x + 12*pow(x,4))*special.erfc(x))/(12.*np.sqrt(2))
    else:
        print("error: index not recognized")
    print("error: harmonic not recognized")