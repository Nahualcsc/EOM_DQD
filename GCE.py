import numpy as np
from inputs import *

class GCE_Impurity(object):
    """Class describing 2-orbital Anderson impurity"""

    def __init__(self,U1,U2,U12):
        self.U1 = U1
        self.U2 = U2
        self.U12 = U12


    def set_gate(self,v1,v2):
        self.v1 = v1
        self.v2 = v2
        
    def energy(self, n1u, n1d, n2u, n2d ):
        v1 = self.v1
        v2 = self.v2 
        n1 = n1u + n1d
        n2 = n2u + n2d
        E = v1*n1 + v2*n2 + self.U1*n1u*n1d + self.U2*n2u*n2d + self.U12*n1*n2 
        return E

    # Compute density (n1,n2) 
    def comp_dens(self,beta):
        Egs=1e+10
        for n1u in range(2):
            for n1d in range(2):
                for n2u in range(2):
                    for n2d in range(2):
                        E=self.energy(n1u,n1d,n2u,n2d)
                        if E<Egs: Egs=E
        
        (n1,n2) = (0.,0.)
        Z = 0.
        for n1u in range(2):
            for n1d in range(2):
                for n2u in range(2):
                    for n2d in range(2):
                        E=self.energy(n1u,n1d,n2u,n2d)-Egs
                        n1 += (n1u+n1d)*np.exp(-beta*E)
                        n2 += (n2u+n2d)*np.exp(-beta*E)
                        Z += np.exp(-beta*E)
        n1 /= Z
        n2 /= Z
        return [n1,n2]
