import numpy as np
from scipy.special import psi
import scipy.optimize as optimize
import scipy.integrate as integrate
from inputs import *

class Equations_of_motion(object):
    def __init__(self, v1,v2,V,TL,TR,U1,U2,U12,gamma1,gamma2):
        self.v1 = v1
        self.v2 = v2
        self.TL = TL
        self.TR = TR
        if TL==TR: self.T = TL
        self.VL =V/2
        self.VR = -V/2
        self.U1 = U1
        self.U2 = U2
        self.U12 = U12
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        
        [self.n1,self.n2] = self.densities()

    ###############################################################
    #### G0 : Diagonal uncontacted system
    ####      Atntti-Pecka like broadening
    ###############################################################
    def int_dig(self,pole,gamma):
        if gamma>0.005 and T>0.005:
            zsL = 0.5+gamma/(4.*np.pi*self.TL)+1.j*(pole-self.VL)/(2.*np.pi*self.TL)
            zsR = 0.5+gamma/(4.*np.pi*self.TR)+1.j*(pole-self.VR)/(2.*np.pi*self.TR)
            integral =  0.5 -(psi(zsL).imag+psi(zsR).imag)/(2.*np.pi)
        elif gamma<0.01 and T>0.001:
             #    #print('EOM with fermi')
            integralL = 1./ (np.exp((pole-self.VL)/self.TL)+1.)
            integralR = 1./ (np.exp((pole-self.VR)/self.TR)+1.)
            integral = 0.5*(integralL+integralR)
        elif gamma>0.01 and T<0.01:
        #    #print('EOM with arctan')
            integral = (1./(np.pi))*(-1.*np.arctan(2*pole/gamma))+0.5
        return integral

    def sys_dens(self,gate, gatebar, U, Ubar, gamma, gammabar):
        fa = self.int_dig(gate              ,gamma)
        fb = self.int_dig(gate + U         ,gamma)
        fc = self.int_dig(gate + U +    self.U12,gamma)
        fd = self.int_dig(gate + U + 2.*self.U12,gamma)
        fe = self.int_dig(gate      +    self.U12,gamma)
        ff = self.int_dig(gate      + 2.*self.U12,gamma)
        fbb =self.int_dig(gatebar + Ubar         ,gammabar)
        fcb =self.int_dig(gatebar + Ubar +    self.U12,gammabar)
        fdb =self.int_dig(gatebar + Ubar + 2.*self.U12,gammabar)
        g = ff/(1.-(fd-ff))
        h = 1.-g*(fbb+fdb-2.*fcb)*(fd-fc)
        A = (fbb*(-2. + fc) - 2.*(-1. + fc)*fcb + fc*fdb)*(fe - ff + (fc - fd - fe + ff)*g) + (1. - fc + fe)*h
        B =(fbb*(-1 + fc)*(-1. + fe) - 2.*(-1. + fc)*fcb*fe + fc*fdb*fe)/A
        C =(fbb*(ff + fe*(-1. + g) - (fc - fd + ff)*g) + fe*h)/A

        tau_alpha = fa
        eta_alphaalpha = 1.+fa-fb
        eta_alphabaralpha = 2.*(fa-fe)+B*(2.*fe-fa-ff)+2.*C*(fb+fe-fa-fc)+g*B*(2.*fc+fa+ff-2.*fe-fb-fd)
        return tau_alpha,eta_alphaalpha,eta_alphabaralpha

    def densities(self):
        tau1, eta11, eta12 = self.sys_dens(self.v1,self.v2,self.U1,self.U2,self.gamma1,self.gamma2)
        tau2, eta22, eta21 = self.sys_dens(self.v2,self.v1,self.U2,self.U1,self.gamma2,self.gamma1)
        denom = 1./(eta11*eta22-eta12*eta21)
        dens1 = denom*(tau1*eta22-tau2*eta12)
        dens2 = denom*(tau2*eta11-tau1*eta21)
        return [dens1,dens2]


    def correlators(self,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar):
        fa = self.int_dig(gate              ,gamma)
        fb = self.int_dig(gate + U         ,gamma)
        fc = self.int_dig(gate + U +    self.U12,gamma)
        fd = self.int_dig(gate + U + 2.*self.U12,gamma)
        fe = self.int_dig(gate      +    self.U12,gamma)
        ff = self.int_dig(gate      + 2.*self.U12,gamma)
        fbb =self.int_dig(gatebar + Ubar         ,gammabar)
        fcb =self.int_dig(gatebar + Ubar +    self.U12,gammabar)
        fdb =self.int_dig(gatebar + Ubar + 2.*self.U12,gammabar)
        g = ff/(1.-(fd-ff))
        h = 1.-g*(fbb+fdb-2.*fcb)*(fd-fc)
        A = (fbb*(-2. + fc) - 2.*(-1. + fc)*fcb + fc*fdb)*(fe - ff + (fc - fd - fe + ff)*g) + (1. - fc + fe)*h
        B =(fbb*(-1 + fc)*(-1. + fe) - 2.*(-1. + fc)*fcb*fe + fc*fdb*fe)/A
        C =(fbb*(ff + fe*(-1. + g) - (fc - fd + ff)*g) + fe*h)/A
        _one_ = nalpha
        _two_ = nalphabar
        _three_ = B*_two_
        _four_ = C*_two_
        _six_ = g*_three_
        _seven_ = (fd-fc)*g*B*_two_ + fc*C*_two_
        return _three_, _four_, _six_, _seven_

    def numerators(self,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar):
        _three_, _four_, _six_, _seven_ = self.correlators(gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar)
        z1 = (1-nalpha-2*nalphabar+_three_+2*_four_-_six_)
        z2 = (nalpha-2*_four_+_six_)
        z3 = (2*_four_-2*_six_)
        z4 = (_six_)
        z5 = (2*nalphabar-2*_three_-2*_four_+2*_six_)
        z6 = (_three_-_six_)
        return [z1,z2,z3,z4,z5,z6]



    def G0(self,w,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar, Delta):
        wa = 1./(w-gate+Delta)
        wb = 1./(w-gate-U+Delta)
        wc = 1./(w-gate-U-self.U12+Delta)#/2
        wd = 1./(w-gate-U-2*self.U12+Delta)
        we = 1./(w-gate-self.U12+Delta)#/2
        wf = 1./(w-gate-2*self.U12+Delta)
        [z1, z2, z3, z4, z5, z6] = self.numerators(gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar)
        return z1*wa+ z2*wb + z3*wc + z4*wd + z5*we + z6*wf


    def Ai(self,w,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar):
        Gr = self.G0( w,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar, 1j*gamma/2)
        return -(Gr).imag/np.pi
    

    def int_dig_I(self,pole,gamma): 
        zsL = 0.5+gamma/(4.*np.pi*self.TL)+1.j*(pole-self.VL)/(2.*np.pi*self.TL)
        zsR = 0.5+gamma/(4.*np.pi*self.TR)+1.j*(pole-self.VR)/(2.*np.pi*self.TR)
        integral = 2*gamma*(-psi(zsL).imag+psi(zsR).imag)/(np.pi) # signs are correct
        return integral

    def int_dig_Q(self,pole,gamma): # over gamma 
        I = self.int_dig_I(pole,gamma)
        zsL = 0.5+gamma/(4.*np.pi*self.TL)+1.j*(pole-self.VL)/(2.*np.pi*self.TL)
        zsR = 0.5+gamma/(4.*np.pi*self.TR)+1.j*(pole-self.VR)/(2.*np.pi*self.TR)
        psi_T = (self.TL-self.TR)/(0.5*(self.TL+self.TR))
        W=(gamma**2/(np.pi))*(psi(zsL).real -psi(zsR).real + np.log(1.+psi_T)-np.log(1.-psi_T))# *gamma
        return (W +(pole-self.VL)*I)



    def I_alpha(self,gate, U,gamma,nalpha,numerators):
        [z1, z2, z3, z4, z5, z6] = numerators
        wa = self.int_dig_I(gate              ,gamma)
        wb = self.int_dig_I(gate + U         ,gamma)
        wc = self.int_dig_I(gate + U  +    self.U12,gamma)
        wd = self.int_dig_I(gate + U  + 2.*self.U12,gamma)
        we = self.int_dig_I(gate      +    self.U12,gamma)
        wf = self.int_dig_I(gate      + 2.*self.U12,gamma)
        return 0.25*(z1*wa + z2*wb + z3*wc + z4*wd + z5*we + z6*wf)

    def Q_alpha(self,gate, U,gamma,nalpha,numerators):
        [z1, z2, z3, z4, z5, z6] = numerators
        wa = self.int_dig_Q(gate              ,gamma)
        wb = self.int_dig_Q(gate + U         ,gamma)
        wc = self.int_dig_Q(gate + U  +    self.U12,gamma)
        wd = self.int_dig_Q(gate + U  + 2.*self.U12,gamma)
        we = self.int_dig_Q(gate      +    self.U12,gamma)
        wf = self.int_dig_Q(gate      + 2.*self.U12,gamma)
        return 0.25*(z1*wa + z2*wb + z3*wc + z4*wd + z5*we + z6*wf)

    def Current(self,n1,n2): # over gamma 
        nums1 = self.numerators(self.v1, self.v2, self.U1, self.U2, self.gamma1, self.gamma2, n1, n2)
        nums2 = self.numerators(self.v2, self.v1, self.U2, self.U1, self.gamma2, self.gamma1, n2, n1)
        I = self.I_alpha(self.v1,self.U1, self.gamma1, n1, nums1)+\
            self.I_alpha(self.v2,self.U2, self.gamma2, n2, nums2)
        Q = self.Q_alpha(self.v1,self.U1, self.gamma1, n1, nums1)+\
            self.Q_alpha(self.v2,self.U2, self.gamma2, n2, nums2)
        return I/self.gamma1,Q/self.gamma1



    ###############################################################
    #### Hartree like decoupling scheme. G(w) and related densities
    ###############################################################

    def fermi (self,pole, T):
        return  1./ (np.exp(pole/T)+1.)

    def G_Hartree(self,w,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar):
        sigma = -0.5j*gamma
        gate = gate+2*U12*nalphabar
        return (1. - nalpha)/(w -gate -sigma*(1.+nalpha*U/(w-gate-U))) + nalpha/(w -gate-U -sigma*(1.-(1.-nalpha)*U/(w-gate))) 


    def tau_Hartree_(self,w,n1,n2):
        G1 = self.G_Hartree( w, self.v1, self.v2, self.U1, self.U2, self.gamma1, self.gamma2,n1,n2)
        G2 = self.G_Hartree( w, self.v2, self.v1, self.U2, self.U1, self.gamma2, self.gamma1,n2,n1)
        return pow(0.5*self.gamma1,2)*(np.conj(G1)*G1).real+pow(0.5*self.gamma2,2)*(np.conj(G2)*G2).real
    def tau_Hartree(self,w,n1,n2):
        A_1 = self.Ai_Hartree(w,self.v1, self.U1, self.gamma1,self.n1,self.n2)
        A_2 = self.Ai_Hartree(w,self.v2, self.U2, self.gamma2, self.n2,self.n1)
        A0 = (A_1+A_2)
        return 0.5*(A0)*np.pi*gamma


    def Ai_Hartree(self,w, gate, U, gamma,nalpha,nalphabar):
        etabroad = 0.0001j
        gate = gate+2*U12*nalphabar
        X=(1.+nalpha*U/(w-gate-U+etabroad))
        Y = (1.-(1.-nalpha)*U/(w-gate+etabroad))
        P1 = -(1-nalpha)*X*0.5*gamma/(pow(w-gate,2)+pow(0.5*gamma*X,2))
        P2 = -nalpha*Y*0.5*gamma/(pow(w-gate-U,2)+pow(0.5*gamma*Y,2))
        Ai= -1/np.pi*(P1+P2)
        return Ai.real

    def n_i_integrand(self,w, gate, U, nalpha,nalphabar):
        integrand = 0.5*(self.fermi(w-self.VL,self.TL)+self.fermi(w-self.VR,self.TR))*self.Ai_Hartree(w, gate, U, gamma,nalpha,nalphabar)
        return integrand

    def solve_H(self, densH): 
        n0, err = integrate.quad(self.n_i_integrand, -np.infty, np.infty, args =( self.v1, U1,densH[0],densH[1]))
        n1, err = integrate.quad(self.n_i_integrand, -np.infty, np.infty, args =( self.v2, U2,densH[1],densH[0]))
        return densH[0] - n0.real, densH[1] - n1.real


    def I_int_EOMH(self,w,n1,n2):
        return 1/(np.pi)*self.tau_Hartree(w,n1,n2)*(self.fermi((w-self.VL),self.TL)-self.fermi((w-self.VR),self.TR))
    def Q_int_EOMH(self,w,n1,n2):
        return 1/(np.pi)*self.tau_Hartree(w,n1,n2)*(w-self.VL)*(self.fermi((w-self.VL),self.TL)-self.fermi((w-self.VR),self.TR))

    def Current_EOMH(self,n1,n2):# over gamma 
        I, err = integrate.quad(self.I_int_EOMH, -np.infty, np.infty,  args =(n1,n2))
        Q, err = integrate.quad(self.Q_int_EOMH, -np.infty, np.infty,  args =(n1,n2))
        return I/gamma,Q/gamma



    def efficiency(self,I,Q):
        output_power=-I*(self.VL-self.VR)
        QR =-Q+output_power
        Q_tot=[Q,QR]
        Qpos = 0.
        for Qi in Q_tot:
            if Qi>0.:
                Qpos +=Qi 
        if Qpos!=0. and output_power>0:
            efficiency = output_power/Qpos
        else: efficiency = 0.
        Q_tot.remove(QR)
        ## CARNOT EFFICIENCY
        QR_carnot=-self.TR*Q/self.TL
        Q_tot.append(QR_carnot)
        Qpos = 0.
        for Qi in Q_tot:
            if Qi>0.:
                Qpos +=Qi 
        if Qpos!=0.:
            c_eff = abs(sum(Q_tot)/Qpos)
        else: c_eff = 1.
        return efficiency/c_eff, output_power


    def L11_integrand_H(self,w,n1,n2):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau_Hartree(w,n1,n2)*fermiprime

    def L12_integrand_H(self,w,n1,n2):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau_Hartree(w,n1,n2)*w*fermiprime

    def L22_integrand_H(self,w,n1,n2):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau_Hartree(w,n1,n2)*pow(w,2)*fermiprime
    

    def transport_coeffs_ints_H(self,n1,n2):
        L11, err = integrate.quad(self.L11_integrand_H, -np.infty, np.infty, args= (n1,n2))
        L12, err = integrate.quad(self.L12_integrand_H, -np.infty, np.infty, args= (n1,n2))
        L22, err = integrate.quad(self.L22_integrand_H, -np.infty, np.infty, args= (n1,n2))
        kappa = L22-pow(L12,2)/L11
        return L11,-(L12/L11)/T,kappa/T

    def tau(self,w):
        A_1 = self.Ai(w,self.v1, self.v2, self.U1, self.U2, self.gamma1, self.gamma2,self.n1,self.n2)
        A_2 = self.Ai(w,self.v2, self.v1, self.U2, self.U1, self.gamma2, self.gamma1,self.n2,self.n1)
        A0 = (A_1+A_2)
        return 0.5*(A0)*np.pi*gamma

    def L11_integrand(self,w):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau(w)*fermiprime

    def L12_integrand(self,w):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau(w)*w*fermiprime

    def L22_integrand(self,w):
        fermiprime = (-1./(2.*self.T))/(np.cosh(w/self.T)+1.)
        return -(1./np.pi)*self.tau(w)*pow(w,2)*fermiprime
    

    def transport_coeffs_ints(self):
        L11, err = integrate.quad(self.L11_integrand, -np.infty, np.infty)
        L12, err = integrate.quad(self.L12_integrand, -np.infty, np.infty)
        L22, err = integrate.quad(self.L22_integrand, -np.infty, np.infty)
        kappa = L22-pow(L12,2)/L11
        return L11,-(L12/L11)/T,kappa/T
    ##################################################################
    ######### NUMERICAL RESOLUTION FOR COUPLED SYSTEM ################
    ######### TO BE COMPLETED ########################################
    ##################################################################

    def G01_coupled_aux(self,w):
        #_five_ dont contribute
        Delta = 0.0001j
        z1 = (1-self._1_-2*self._2_+self._22_+2*self._12_-self._122_)
        z2 = (self._1_-2*self._12_+self._122_)
        z3 = (2*self._12_-2*self._122_)
        z4 = self._122_
        z5 = (2*self._2_-2*self._22_-2*self._12_+2*self._122_)
        z6 = (self._22_-self._122_)
        wa = 1./(w-self.v1+Delta)
        wb = 1./(w-self.v1-self.U1+Delta)
        wc = 1./(w-self.v1-self.U1-self.U12+Delta)
        wd = 1./(w-self.v1-self.U1-2*self.U12+Delta)
        we = 1./(w-self.v1-self.U12+Delta)
        wf = 1./(w-self.v1-2*self.U12+Delta)
        G0 = z1*wa+ z2*wb + z3*wc + z4*wd + z5*we + z6*wf
        return G0.real
    def G02_coupled_aux(self,w):
        #_five_ dont contribute
        Delta = 0.0001j
        z1 = (1-self._2_-2*self._1_+self._11_+2*self._12_-self._211_)
        z2 = (self._2_-2*self._12_+self._211_)
        z3 = (2*self._12_-2*self._211_)
        z4 = self._211_
        z5 = (2*self._1_-2*self._11_-2*self._12_+2*self._211_)
        z6 = (self._11_-self._211_)
        wa = 1./(w-self.v2+Delta)
        wb = 1./(w-self.v2-self.U2+Delta)
        wc = 1./(w-self.v2-self.U2-self.U12+Delta)
        wd = 1./(w-self.v2-self.U2-2*self.U12+Delta)
        we = 1./(w-self.v2-self.U12+Delta)
        wf = 1./(w-self.v2-2*self.U12+Delta)
        G0 = z1*wa+ z2*wb + z3*wc + z4*wd + z5*we + z6*wf
        return G0.real
    def fermieff(self,w):
        return 0.5*(self.fermi((w-self.VL),self.TL)+self.fermi((w-self.VR),self.TR))

    def _1approx_int(self,w,n1,n2):
        G0 = self.G0( w,self.v1, self.v2, self.U1, self.U2, self.gamma1, self.gamma2,n1,n2, 0.0001j).real
        Gr = ((G0)/(1+G0*0.5j*gamma))
        A1 = -(Gr).imag/np.pi
        #G_imag =((G01)/(1+G01*0.5j*gamma)).imag
        return self.fermieff(w)*A1

    def _2approx_int(self,w,n1,n2):
        G0 = self.G0( w,self.v2, self.v1, self.U2, self.U1, self.gamma2, self.gamma1,n2,n1, 0.0001j).real
        Gr = ((G0)/(1+G0*0.5j*gamma))
        A2 = -(Gr).imag/np.pi
        #G_imag = -pow(G02*0.5*gamma,2)/(1+pow(G02*0.5*gamma,2))
        #G_imag =((G02)/(1+G02*0.5j*gamma)).imag
        return self.fermieff(w)*A2

    def _1_int(self,w):
        G01 = self.G01_coupled_aux(w)
        G_imag =((G01)/(1+G01*0.5j*gamma)).imag
        return -1./(np.pi)*self.fermieff(w)*G_imag

    def _2_int(self,w):
        G02 = self.G02_coupled_aux(w)
        G_imag = ((G02)/(1+G02*0.5j*gamma)).imag
        return -1./(np.pi)*self.fermieff(w)*G_imag

    def _11_int(self,w):
        G02 = self.G02_coupled_aux(w)
        G_imag = (1./(1+G02*0.5j*gamma)).imag
        wpb2 = w-(self.v1+self.U1)
        wpb3 = w-(self.v1+self.U1+self.U12)
        wpb4 = w-(self.v1+self.U1+2*self.U12)
        factor = (self._1_*wpb2+2*(wpb3-wpb2)*self._12_+(wpb2+wpb4-2*wpb3)* self._2s2sb1s_)
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def _22_int(self,w):
        G01 = self.G01_coupled_aux(w)
        G_imag =((G01)/(1+G01*0.5j*gamma)).imag
        wpb2 = w-(self.v2+self.U2)
        wpb3 = w-(self.v2+self.U2+self.U12)
        wpb4 = w-(self.v2+self.U2+2*self.U12)
        factor = (self._2_*wpb2+2*(wpb3-wpb2)*self._12_+(wpb2+wpb4-2*wpb3)* self._1s1sb2sb_)
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def _12_int(self,w):
        G02 = self.G02_coupled_aux(w)
        G_imag = (1./(1+G02*0.5j*gamma)).imag
        wp3 = w-(self.v1+self.U1+self.U12)
        wp4 = w-(self.v1+self.U1+2*self.U12)
        wp5 = w-(self.v1+self.U12)
        wp6 = w-(self.v1+2*self.U12)
        factor = self._2_*wp5+(wp6-wp5)*self._22_+(wp3-wp5)* self._12_+(wp4+wp5-wp3-wp6)*self._122_
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def _122_int(self,w):
        G01 = self.G01_coupled_aux(w)
        G_imag =((G01)/(1+G01*0.5j*gamma)).imag
        wp4 = w-(self.v1+self.U1+2*self.U12)
        wp6 = w-(self.v1+2*self.U12)
        factor = ((wp4-wp6)*self._122_+(wp6)*self. _22_)
        return -1./(4*np.pi)*self.fermieff(w)*G_imag*factor

    def _211_int(self,w):
        G02 = self.G02_coupled_aux(w)
        G_imag = (1./(1+G02*0.5j*gamma)).imag
        wp4 = w-(self.v2+self.U2+2*self.U12)
        wp6 = w-(self.v2+2*self.U12)
        factor = ((wp4-wp6)*self._211_+(wp6)* self._11_)
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def _2s2sb1s_int(self,w):
        G02 = self.G02_coupled_aux(w)
        G_imag = (1./(1+G02*0.5j*gamma)).imag
        wp4 = w-(self.v2+self.U2+2*self.U12)
        wp3 = w-(self.v2+self.U2+self.U12)
        factor = ((wp4-wp3)*self._211_+(wp3)* self._12_)
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def _1s1sb2sb_int(self,w):
        G01 = self.G01_coupled_aux(w)
        G_imag =((G01)/(1+G01*0.5j*gamma)).imag
        wp4 = w-(self.v1+self.U1+2*self.U12)
        wp3 = w-(self.v1+self.U1+self.U12)
        factor = ((wp4-wp3)*self._122_+(wp3)* self._12_)
        return -1./(np.pi)*self.fermieff(w)*G_imag*factor

    def Solve_coupled_system(self, corr):
        self.Delta = -0.5j*gamma
        self._1_ = corr[0] # _1_ = <n1> 
        self._2_ = corr[1] # _2_ = <n2> 
        self._11_ = corr[2] # _11_ = <n1n1> 
        self._22_ = corr[3] # _22_ = <n2n2> 
        self._12_ = corr[4] # _12_ = <n1n2> 
        self._122_ = corr[5] # _122_ = <n1n2n2> 
        self._211_ = corr[6] # _211_ = <n2n1n1> 
        self._2s2sb1s_ = corr[7] # _221_ = <n2n2n1> 
        self._1s1sb2sb_ = corr[8] # _112_ = <n1n1n2> 
        __1__ ,err = integrate.quad(self._1_int, -np.infty, np.infty)
        __2__ ,err = integrate.quad(self._2_int, -np.infty, np.infty)
        __11__ ,err = integrate.quad(self._11_int, -np.infty, np.infty)
        __22__ ,err = integrate.quad(self._22_int, -np.infty, np.infty)
        __12__ ,err = integrate.quad(self._12_int, -np.infty, np.infty)
        __122__ ,err = integrate.quad(self._122_int, -np.infty, np.infty)
        __211__ ,err = integrate.quad(self._211_int, -np.infty, np.infty)
        __2s2sb1s__ ,err = integrate.quad(self._2s2sb1s_int, -np.infty, np.infty)
        __1s1sb2sb__ ,err = integrate.quad(self._1s1sb2sb_int, -np.infty, np.infty)
        return __1__-corr[0],__2__-corr[1],__11__-corr[2],__22__-corr[3],__12__-corr[4],__122__-corr[5]\
             ,__211__-corr[6],__2s2sb1s__-corr[7],__1s1sb2sb__-corr[8]

    def dens_approx_G(self, corr):
        n1 ,err = integrate.quad(self._1approx_int, -np.infty, np.infty, args = (corr[0],corr[1]), epsabs=1.49e-10, epsrel=1.49e-10)
        n2 ,err = integrate.quad(self._2approx_int, -np.infty, np.infty, args = (corr[0],corr[1]), epsabs=1.49e-10, epsrel=1.49e-10)
        return n1 - corr[0], n2 - corr[1]



    def solve_vs(self,vs,n1_int,n2_int):
        tau1, eta11, eta12 = self.sys_dens(vs[0],vs[1],0.,0.,self.gamma1,self.gamma2)
        tau2, eta22, eta21 = self.sys_dens(vs[1],vs[0],0.,0.,self.gamma2,self.gamma1)
        denom = 1./(eta11*eta22-eta12*eta21)
        dens1 = denom*(tau1*eta22-tau2*eta12)
        dens2 = denom*(tau2*eta11-tau1*eta21)
        return [dens1-n1_int,dens2-n2_int]



    ##################################################################
    ################### EXTRA. FOR DEBUGGING  ########################
    ##################################################################
'''
    #### checked, it gives the same as self.numerators
    def numerators2(self,gate, gatebar, U, Ubar, gamma, gammabar,nalpha,nalphabar):
        fa = self.int_dig(gate              ,gamma)
        fb = self.int_dig(gate + U         ,gamma)
        fc = self.int_dig(gate + U +    self.U12,gamma)
        fd = self.int_dig(gate + U + 2.*self.U12,gamma)
        fe = self.int_dig(gate      +    self.U12,gamma)
        ff = self.int_dig(gate      + 2.*self.U12,gamma)
        fbb =self.int_dig(gatebar + Ubar         ,gammabar)
        fcb =self.int_dig(gatebar + Ubar +    self.U12,gammabar)
        fdb =self.int_dig(gatebar + Ubar + 2.*self.U12,gammabar)
        g = ff/(1.-(fd-ff))
        h = 1.-g*(fbb+fdb-2.*fcb)*(fd-fc)
        A = (fbb*(-2. + fc) - 2.*(-1. + fc)*fcb + fc*fdb)*(fe - ff + (fc - fd - fe + ff)*g) + (1. - fc + fe)*h
        B =(fbb*(-1 + fc)*(-1. + fe) - 2.*(-1. + fc)*fcb*fe + fc*fdb*fe)/A
        C =(fbb*(ff + fe*(-1. + g) - (fc - fd + ff)*g) + fe*h)/A
        z1 = 1-nalpha+nalphabar*(B+2*C-g*B-2)
        z2 = nalpha+nalphabar*(g*B-2*C)
        z3 = 2*nalphabar*(C-g*B)
        z4 = nalphabar*(g*B)
        z5 = 2*nalphabar*(1-B-C+g*B)
        z6 = nalphabar*(1-g)*B
        return [z1,z2,z3,z4,z5,z6]
'''

    ###### numerical implementation of the charge and heat currents
    ###### gives the same as the analytical implementation
    ###### I commented out this piece since its slower
'''
    def tau(self,w,n1,n2):
        A_1 = self.Ai(w,self.v1, self.v2, self.U1, self.U2, self.gamma1, self.gamma2,n1,n2)
        A_2 = self.Ai(w,self.v2, self.v1, self.U2, self.U1, self.gamma2, self.gamma1,n2,n1)
        A0 = (A_1+A_2)
        return (A0)*np.pi*gamma
    def I_integrand(self,w,n1,n2):
        return (1./np.pi)*self.tau(w,n1,n2)*(self.fermi((w-self.VL),self.TL)-self.fermi((w-self.VR),self.TR))
    def Q_integrand(self,w,n1,n2):
        return (1./np.pi)*(w-self.VL)*self.tau(w,n1,n2)*(self.fermi((w-self.VL),self.TL)-self.fermi((w-self.VR),self.TR))
    def Current_int(self,n1,n2):
        I, err = integrate.quad(self.I_integrand, -np.infty, np.infty,  args =(n1,n2))
        Q, err = integrate.quad(self.Q_integrand, -np.infty, np.infty,  args =(n1,n2))
        return I/gamma,Q/gamma
'''
