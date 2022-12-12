# -*- coding: utf-8 -*-
"""
Computes H functions, for any integer order n.
Not efficient by any means, just for testing.
You can use e.g. Boost's hypergeometric pFq in C++ for more performant evaluation.
"""

import mpmath as mp

mp.mp.dps = 100  # set precision, probably overkill for the H functions.

#helper functions, used to define the actual H functions, see paper Machielsen
def helper1(x, n):
    if n == 0:
        return 0  # edge case, value irrelevant because it will be multiplied by 0 anyway.
    else:
        n = abs(n)
        return mp.hyp2f2(-n, n, 1/2, 1, -x*x)/(2*n)-2*x*mp.hyp2f2(1/2-n, 1/2+n, 3/2, 3/2, -x*x)/mp.sqrt(mp.pi)


def helper2(x, n):
    n = abs(n)
    return mp.hyp2f2(1/2-n, 1/2+n, 1/2, 3/2, -x*x)/mp.sqrt(2*mp.pi) - 1/2*mp.sqrt(2)*n*x*mp.hyp2f2(1-n, 1+n, 3/2, 2, -x*x)


def helper3(x, n):
    n = abs(n)
    return mp.hyp2f2(1/2-n, 1/2+n, 1/2, 1/2, -x*x)/(4*mp.sqrt(mp.pi)*x) - 1/2*n*mp.hyp2f2(1-n, 1+n, 1, 3/2, -x*x)


def helper4(x, n):
    n = abs(n)
    return n*(n*n-1)*x/(3*mp.sqrt(2))*mp.hyp2f2(2-n, 2+n, 2, 5/2, -x*x) + mp.hyp2f2(1/2-n, 1/2+n, -1/2, 1/2, -x*x)/(8*mp.sqrt(2*mp.pi)*x*x)


def helper5(x, n):
    n = abs(n)
    return n*(n*n-1)/6*mp.hyp2f2(2-n, 2+n, 1, 5/2, -x*x) - mp.hyp2f2(1/2-n, 1/2+n, -1/2, -1/2, -x*x)/(32*mp.sqrt(mp.pi)*x*x*x)


def helper6(x, n):
    n = abs(n)
    return 3/(32*mp.sqrt(mp.pi)*x*x*x)*mp.hyp2f2(1/2-n, 1/2+n, -3/2, 1/2, -x*x) - n*(4-5*n*n+n**4)/30*x*x*mp.hyp2f2(3-n, 3+n, 3, 7/2, -x*x)


def getH(x, index, n, jacobian=False):
    '''
    Parameters
    ----------
    x : radial position, so >= 0
    index : select H0, H0breve, H1, H2, H3, H4, H5
    n : harmonic, integer
    jacobian : bool, optional
        if True, multiply by x so that the singularity is avoided.
        The default is False.

    Returns
    -------
    double
    '''
    if x == 0:
        return getH(1.0e-15, index, n, jacobian) #Evaluation at exactly x=0 might result in issues.

    ans = 0.0
    if index == 0:
        ans = helper5(x, n) - (helper5(x, n-1)+helper5(x, n+1))/2 #H0
    elif index == 1:
        ans = helper6(x, n) - (helper6(x, n-1)+helper6(x, n+1))/2 #H0breve
    elif index==2:
        ans = n*n*helper1(x, n) #H1
    elif index==3:
        ans = n*(helper3(x, n) - (helper3(x, n-1)+helper3(x, n+1))/2) #H2
    elif index==4:
        ans = helper3(x,n) #H3
    elif index==5:
        ans = n*helper2(x,n) #H4
    else:
        ans = helper4(x, n) - (helper4(x, n-1)+helper4(x, n+1))/2 #H5

    if jacobian:
        return ans*x
    else:
        return ans
