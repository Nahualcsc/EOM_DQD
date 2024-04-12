import numpy as np


# bummer -- no complex trigamma in scipy.special!
# roll our own from a continued fraction

def trigamma(z,eps=np.spacing(np.double(1)),realmin=np.finfo(np.double).tiny):
    """ Trigramma function of a complex argument for arg(z) < pi/2.
The trigramma and tetragamma functions are the only two polygamma functions
for which closed-form numerator coeffients are known [2].

References:
[1] Cuyt et al., _Handbook of Continued Fractions for Special Functions_.
[2] Catherine M. Bonan-Hamada and William B. Jones.
"Stieltjes continued fractions for polygamma functions; speed of convergence,"
J. Computational and Appl. Math. 179(1--2):47--55 (2005)
<http://www.sciencedirect.com/science/article/pii/S037704270400442X>.
"""
    # Use trigamma reflection on the left-half plane
    zcomplex = isinstance(z,(complex,np.complex128,np.complex64))
    if (zcomplex and z.real < 0) or (not zcomplex and z < 0):
        return -trigamma(1-z) + (np.pi/np.sin(np.pi*z))**2
    if z == 0: return np.inf
    # Use basic trigamma recurrence for small arguments
    if np.abs(z) < 10.5: return trigamma(z+1) + 1/(z*z)
    zinv = 1/z
    return zinv*(1 + 0.5*zinv + 2*np.pi*trigamma_g1(z,eps=eps,realmin=realmin))

def trigamma_g1(z,maxits=40,eps=np.spacing(np.double(1)),realmin=np.finfo(np.double).tiny):
    z2 = z*z
    a1 = 1/(12*np.pi)
    def an(n):
        """ Defined for n > 1. """
        n2 = np.double(n)*n
        return n2/(4*n2-1)*(n2-1)/4
    # adaptation of NR's Lentz's continued fraction algorithm
    # http://www.aip.de/groups/soe/local/numres/bookfpdf/f5-2.pdf
    # hardcode n = 1 initialization
    a = a1
    b = z2
    d = 1/b
    if np.abs(d) < realmin: d = realmin
    c = a/realmin
    cf = a/b
    for n in range(2,maxits+1):
        a = an(n)
        b = 1 if n%2 == 0 else z2
        d = a*d + b
        if np.abs(d) < realmin: d = realmin
        c = b + a/c
        if np.abs(c) < realmin: c = realmin
        d = 1/d
        dlt = d*c
        cf *= dlt
        if np.abs(dlt-1) < eps: break
    else:
        raise Exception('trigamma_g1 failed to converge.')
    return cf
