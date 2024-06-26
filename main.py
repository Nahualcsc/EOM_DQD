import numpy as np
from scipy import integrate
from tqdm import tqdm
import os

from inputs import *
from EOM import *
from GCE import *
from generate_plots import *

#######################################################
#######################################################
Uave = 0.5*(U1+U2)

try:
    if U12<=U1 and U12<=U2: ##REGIME 1
        print('REGIME 1')
    elif U12>U1 and Uave>U12: #REGIME 2
        print('REGIME 2')
    elif Uave>=U1 and U12>=Uave: #REGIME 3
        print('REGIME 3')
except: pass
list_file = []

if calculate_densities_vs_gate:
    imp = GCE_Impurity(U1,U2,U12)
    dens_EOMH = [1.,1.]
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g_{}'.format(U1, U2, U12, T, gamma)
    for veq in tqdm(v_range):
        imp.set_gate( veq+delta_v/2,veq-delta_v/2)
        dens_GCE = imp.comp_dens(beta)
        EOM = Equations_of_motion( veq+delta_v/2,veq-delta_v/2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
        I_EOM0,Q_EOM0 = EOM.Current(EOM.n1,EOM.n2)
        if delta_T!=0.:efficiency_EOM0, output_power_EOM0 = EOM.efficiency(I_EOM0,Q_EOM0)
        else: efficiency_EOM0, output_power_EOM0 = 0.,0.
        if compute_EOMH:
            #dens_EOMH = optimize.root( EOM.solve_H, [0.5*dens_GCE[0],0.5*dens_GCE[1]], method='hybr').x
            dens_EOMH = optimize.root( EOM.solve_H, [dens_EOMH[0],dens_EOMH[1]], method='hybr').x
            ##dens_EOMH = optimize.root( EOM.dens_approx_G, [EOM.n1,EOM.n2], method='hybr', options={'maxiter': 500, 'fatol': 1e-12}).x
            ##dens_EOMH = optimize.root( EOM.Solve_coupled_system, [EOM.n1,EOM.n2,0,0,0,0,0,0,0], method='hybr').x
            I_EOMH,Q_EOMH = EOM.Current_EOMH(dens_EOMH[0],dens_EOMH[1])
            if delta_T!=0.:efficiency_EOMH, output_power_EOMH = EOM.efficiency(I_EOMH,Q_EOMH)
            else: efficiency_EOMH, output_power_EOMH,c_eff_EOMH = 0.,0.,0.
            list_file.append((veq,dens_GCE[0],dens_GCE[1],2*dens_EOMH[0],2*dens_EOMH[1],2*EOM.n1,2*EOM.n2,\
                I_EOMH,I_EOM0,Q_EOMH,Q_EOM0,efficiency_EOMH, output_power_EOMH,efficiency_EOM0, output_power_EOM0))
        else: list_file.append((veq,dens_GCE[0],dens_GCE[1],2*EOM.n1,2*EOM.n2, I_EOM0,Q_EOM0,efficiency_EOM0, output_power_EOM0))
    plot_calculate_densities(list_file,output)
    

###############################################
###############################################

if calculate_spectral_function:
    EOM = Equations_of_motion(v1,v2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
    n1 = EOM.n1
    n2 = EOM.n2
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g_{}'.format(U1, U2, U12, T, gamma)
    if compute_EOMH:
        dens_EOMH = optimize.root( EOM.solve_H, [0.5,0.5], method='hybr').x
    for w in w_range: 
        A_1 = EOM.Ai(w,v1, v2, U1, U2, gamma1, gamma2,EOM.n1,EOM.n2)
        A_2 = EOM.Ai(w,v2, v1, U2, U1, gamma2, gamma1,EOM.n1,EOM.n2)
        A0 = (A_1+A_2)
        if compute_EOMH:
            #AH_1 = EOM.Ai_Hartree(w, v1, U1, gamma1,dens_EOMH[0],dens_EOMH[1])
            #AH_2 = EOM.Ai_Hartree(w, v2, U2, gamma2,dens_EOMH[1],dens_EOMH[0])
            #AH = (AH_1+AH_2)
            G01 = EOM.G0( w,v1, v2, U1, U2, gamma1, gamma2,EOM.n1,EOM.n2, 0.0001j).real
            G02 = EOM.G0( w,v2, v1, U2, U1, gamma2, gamma1,EOM.n2,EOM.n1, 0.0001j).real
            AH_1 = -1/(np.pi)*G01/(1+G01*0.5j*gamma).imag
            AH_2 = -1/(np.pi)*G02/(1+G02*0.5j*gamma).imag
            print(AH_2)
            AH = (AH_1+AH_2)
            list_file.append((w,A_1,A_2,A0,AH_1,AH_2,AH))
        else: list_file.append((w,A_1,A_2,A0))
    norma_A1_O, err = integrate.quad(EOM.Ai, -np.infty, np.infty,  args =(v1,v2, U1,U2, gamma1,gamma2,EOM.n1,EOM.n2))
    if compute_EOMH:
        norma_A1_H, err = integrate.quad(EOM.Ai_Hartree, -np.infty, np.infty,  args =(v1, U1, gamma1,dens_EOMH[0],dens_EOMH[1]))
        print('|A1|_H = {}, |A1|_0 = {}'.format(norma_A1_H,norma_A1_O))
    else: print('|A1|_0 = {}'.format(norma_A1_O))
    plot_calculate_spectral_function(list_file,output)



###############################################
###############################################

def compute_Aw(w_range,U12_range,V,delta_v):
    A = np.zeros((len(U12_range), len(w_range)))
    output = 'U1_{}_U2_{}_V_{}_T_{}_g1_{}_g2_{}_dT_{}'.format(U1, U2, V, T, gamma1, gamma2,delta_T)
    for i, U12 in enumerate(U12_range):
        if colormap_spectral_function or movie_spectral_function_vary_V:
            v1 = -0.5*U1-U12
            v2 = -0.5*U2-U12
        elif movie_spectral_function_vary_dv:
            print(delta_v)
            v1 = -0.5*U1-U12+0.5*delta_v
            v2 = -0.5*U2-U12-0.5*delta_v
        for j, w in enumerate(w_range):
            EOM = Equations_of_motion(v1,v2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
            A_1 = EOM.Ai(w,v1, v2, U1, U2, gamma1, gamma2,EOM.n1,EOM.n2)
            A_2 = EOM.Ai(w,v2, v1, U2, U1, gamma2, gamma1,EOM.n2,EOM.n1)
            A[i, j] = (A_1+A_2)
    return A

if colormap_spectral_function:
    output = 'U1_{}_U2_{}_T_{}_g1_{}_g2_{}_dT_{}'.format(U1, U2, T, gamma1, gamma2,delta_T)
    A = compute_Aw(w_range,U12_range,V,delta_v=0)
    plot_colormap_spectral_function(w_range, U12_range, A,'','Aw',output)


if movie_spectral_function_vary_V:
    folder = 'movies'
    output = 'U1_{}_U2_{}_dV_{}_T_{}_g1_{}_g2_{}_dT_{}_V_{}_to_{}'.format(U1, U2, delta_v, T, gamma1, gamma2,delta_T,V_range[0],V_range[-1])
    for V in tqdm(V_range):
        A = compute_Aw(w_range,U12_range,V,delta_v)
        title = '$V$ = '+str(round(V,4))
        plot_colormap_spectral_function(w_range, U12_range, A, title,folder,str(V))
    create_movies_colormap_spectral_function(output)

if movie_spectral_function_vary_dv:
    folder = 'movies'
    output = 'U1_{}_U2_{}_V_{}_T_{}_g1_{}_g2_{}_dT_{}_dv_{}_to_{}'.format(U1, U2, V, T, gamma1, gamma2,delta_T,dv_range[0],dv_range[-1])
    for delta_v in tqdm(dv_range):
        #v1 = v+0.5*delta_v
        #v2 = v -0.5*delta_v
        A = compute_Aw(w_range,U12_range,V,delta_v)
        title = '$\\delta v$ = '+str(round(delta_v,4))
        plot_colormap_spectral_function(w_range, U12_range, A, title,folder,str(delta_v))
    create_movies_colormap_spectral_function(output)

###############################################
###############################################
def compute_colormaps_currents(v_range,V_range,U12):
    dens_0 =  np.zeros((len(V_range), len(v_range)))
    dens_1 =  np.zeros((len(V_range), len(v_range)))
    I_EOM0 =  np.zeros((len(V_range), len(v_range)))
    Q_EOM0 =  np.zeros((len(V_range), len(v_range)))
    for i, V in enumerate(V_range):
        for j, v in enumerate(v_range):
            EOM = Equations_of_motion(v+delta_v/2,v-delta_v/2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            dens = EOM.densities()
            dens_0[i, j] = 2*EOM.n1
            dens_1[i, j] = 2*EOM.n2
            I, Q= EOM.Current(EOM.n1,EOM.n2)
            I_EOM0[i, j] = I
            Q_EOM0[i, j] = Q
    return dens_0,dens_1,I_EOM0,Q_EOM0

if colormaps_currents:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dT_{}_dv_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_T,dv)
    dens_0,dens_1,I_EOM0,Q_EOM0 = compute_colormaps_currents(v_range,V_range,U12)
    plot_colormaps_currents(v_range,V_range,dens_0,dens_1,I_EOM0,Q_EOM0,'','colormaps',output)



if movie_currents_vary_U12:
    output = 'U1_{}_U2_{}_T_{}_DT_{}_g1_{}_g2_{}_dv_{}_U12_{}_to_{}'.format(U1, U2, T, delta_T, gamma1, gamma2,delta_v,U12_range[0],U12_range[-1])
    folder = 'movies'
    for U12 in tqdm(U12_range):
        dens_0,dens_1,I_EOM0,Q_EOM0 = compute_colormaps_currents(v_range,V_range,U12)
        title = '$U_{12}$ = '+str(round(U12,3))
        plot_colormaps_currents(v_range,V_range,dens_0,dens_1,I_EOM0,Q_EOM0,title,folder,str(U12))
    create_movies_colormap_currents(output)


###############################################
###############################################

def compute_colormap_efficiency(v_range,V_range,delta_T):
    eficciency_norm_EOM0 = np.zeros((len(V_range), len(v_range)))  # Initialize the efficiency array
    for i, V in enumerate(V_range):
        for j, v in enumerate(v_range):
            v1 = v + delta_v
            v2 = v - delta_v
            EOM = Equations_of_motion(v1, v2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            I, Q = EOM.Current(EOM.n1, EOM.n2)
            eficciency_norm_EOM0[i, j], _ = EOM.efficiency(I, Q)
    return eficciency_norm_EOM0

if colormap_efficiency:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dT_{}_dv_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_T,dv)
    eficciency_norm_EOM0 = compute_colormap_efficiency(v_range,V_range,delta_T)
    plot_colormap_efficiency(v_range,V_range,eficciency_norm_EOM0,V,U12,'','colormaps',output)


if movie_efficiency_vary_DT:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dv_{}_DT_{}_to_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_v,DT_range[0],DT_range[-1])
    folder = 'movies'
    for delta_T in tqdm(DT_range):
        eficciency_norm_EOM0 = compute_colormap_efficiency(v_range,V_range,delta_T)
        title = '$\\Delta T$ = '+str(round(delta_T,4))
        plot_colormap_efficiency(v_range,V_range,eficciency_norm_EOM0,delta_T,title,folder,'efficiency_'+str(delta_T))
    create_movies_colormap_efficiency(output)

###############################################
###############################################

def compute_SD(v1_range,v2_range,V,U12,delta_T):
    A = np.zeros((len(v1_range), len(v2_range), 6))  # For z1 to z6
    B = np.zeros((len(v1_range), len(v2_range), 6))  # For z1 to z6
    Z = np.zeros((len(v1_range), len(v2_range), 6))  # For z1 to z6
    I_values = np.zeros((len(v1_range), len(v2_range)))
    Q_values = np.zeros((len(v1_range), len(v2_range)))
    nsum = np.zeros((len(v1_range), len(v2_range)))
    for i, v1 in enumerate(v1_range):
        for j, v2 in enumerate(v2_range):
            EOM = Equations_of_motion(v1, v2, V, T + 0.5*delta_T, T -0.5*delta_T, U1, U2, U12,gamma1, gamma2)
            nsum[i,j] = (2*EOM.n1+5*EOM.n2)
            A[i, j, :] = EOM.numerators(v2, v1, U2, U1, gamma2, gamma1, EOM.n2, EOM.n1)
            EOM = Equations_of_motion(v1, v2, V, T + 0.5*delta_T, T -0.5*delta_T, U1, U2, U12, 0.00001, 0.00001)
            B[i, j, :] = EOM.numerators(v2, v1, U2, U1, 0.00001, 0.00001, EOM.n2, EOM.n1)#EOM.numerators(v1, v2, U1, U2, gamma1, gamma2, EOM.n1, EOM.n2)
            Z[i, j, :] = A[i, j, :] -B[i, j, :] 
            I, Q = EOM.Current(EOM.n1, EOM.n2)
            I_values[i, j] = I*100
            Q_values[i, j] = Q*100
    return I_values, Q_values,Z, nsum


if colormap_SD:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dT_{}_V_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_T,V)
    I_values, Q_values,Z, nsum = compute_SD(v1_range,v2_range,V,U12,delta_T)
    plot_colormap_SD(v1_range, v2_range, I_values, Q_values,Z, nsum,V,U12,delta_T,'','colormaps',output)


if movie_SD_vary_V:
    folder = 'movies'
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dT_{}_V_{}_to_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_T,V_range[0],V_range[-1])
    for V in tqdm(V_range):
        I_values, Q_values,Z, nsum = compute_SD(v1_range,v2_range,V,U12,delta_T)
        title = 'V = '+str(round(V,4))
        plot_colormap_SD(v1_range, v2_range, I_values, Q_values,Z, nsum,V,U12,delta_T,title,folder,str(V))
    create_movies_colormap_SD(output)

if movie_SD_vary_U12:
    folder = 'movies'
    output = 'U1_{}_U2_{}_T_{}_g1_{}_g2_{}_V_{}_dT_{}_U12_{}_to_{}'.format(U1, U2, T, gamma1, gamma2,V,delta_T,U12_range[0],U12_range[-1])
    for U12 in tqdm(U12_range):
        I_values, Q_values,Z, nsum = compute_SD(v1_range,v2_range,V,U12,delta_T)
        title = '$U12$ = '+str(round(U12,2))
        plot_colormap_SD(v1_range, v2_range, I_values, Q_values,Z, nsum,V,U12,delta_T,title,folder,str(U12))
    create_movies_colormap_SD(output)

if movie_SD_vary_DT:
    folder = 'movies'
    output = 'U1_{}_U2_{}_T_{}_g1_{}_g2_{}_V_{}_U12_{}_dT_{}_to_{}'.format(U1, U2, T, gamma1, gamma2,V,U12,DT_range[0],DT_range[-1])
    for delta_T in tqdm(DT_range):
        I_values, Q_values,Z, nsum = compute_SD(v1_range,v2_range,V,U12,delta_T)
        title = '$\\Delta T$ = '+str(round(delta_T,4))
        plot_colormap_SD(v1_range, v2_range, I_values, Q_values,Z, nsum,V,U12,delta_T,title,folder,str(delta_T))
    create_movies_colormap_SD(output)

###############################################
###############################################


if calculate_transport_coeffs:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g_{}'.format(U1, U2,U12,T, gamma)
    for veq in tqdm(v_range):
        EOM = Equations_of_motion( veq+delta_v/2,veq-delta_v/2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
        G,S,kappa = EOM.transport_coeffs()
        ZT = T*G*pow(S,2)/kappa
        if compute_EOMH:
            dens_EOMH = optimize.root( EOM.solve_H, [0.5,0.5], method='hybr').x
            G_H,S_H,kappa_H = EOM.transport_coeffs_ints_H(dens_EOMH[0],dens_EOMH[1])
            ZT_H = T*G_H*pow(S_H,2)/kappa_H
            list_file.append((veq,G*np.pi,S,kappa*np.pi,ZT,G_H*np.pi,S_H,kappa_H*np.pi,ZT_H))
        else:
            list_file.append((veq,G*np.pi,S,kappa*np.pi,ZT))
    plot_calculate_transport_coeffs(list_file,output)



def compute_colormaps_transport_coeffs_vary_U12(v_range,U12_range):
    G =  np.zeros((len(U12_range), len(v_range)))
    S =  np.zeros((len(U12_range), len(v_range)))
    kappa =  np.zeros((len(U12_range), len(v_range)))
    ZT =  np.zeros((len(U12_range), len(v_range)))
    for i, U12 in enumerate(tqdm(U12_range)):
        for j, v in enumerate(v_range):
            EOM = Equations_of_motion(v+delta_v/2,v-delta_v/2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            G_,S_,kappa_ = EOM.transport_coeffs()
            ZT_ = T*G_*pow(S_,2)/kappa_
            G[i, j] = G_
            S[i, j] = S_
            kappa[i, j] = kappa_
            ZT[i, j] = ZT_
    return G,S,kappa ,ZT

def compute_colormaps_transport_coeffs_vary_dv(v_range,dV_range):
    G =  np.zeros((len(dv_range), len(v_range)))
    S =  np.zeros((len(dv_range), len(v_range)))
    kappa =  np.zeros((len(dv_range), len(v_range)))
    ZT =  np.zeros((len(dv_range), len(v_range)))
    for i, delta_v in enumerate(tqdm(dv_range)):
        for j, v in enumerate(v_range):
            EOM = Equations_of_motion(v+delta_v/2,v-delta_v/2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            G_,S_,kappa_ = EOM.transport_coeffs()
            ZT_ = T*G_*pow(S_,2)/kappa_
            G[i, j] = G_
            S[i, j] = S_
            kappa[i, j] = kappa_
            ZT[i, j] = ZT_
    return G,S,kappa ,ZT

if colormaps_transport_coeffs_vary_U12:
    output = 'U1_{}_U2_{}_T_{}_g_{}_dv_{}_U12_from_{}_to_{}'.format(U1, U2, T, gamma1,delta_v,U12_range[0],U12_range[-1])
    title = '$U_1=$'+str(U1)+', $U_2=$'+str(U2)+', $\\delta v=$'+str(delta_v)+', $T=$'+str(T)+', $\\gamma=$'+str(gamma)
    G,S,kappa ,ZT = compute_colormaps_transport_coeffs_vary_U12(v_range,U12_range)
    plot_colormaps_transport_coeffs(v_range,U12_range,G,S,kappa ,ZT,title,'$U_{12}$','colormaps',output)

if colormaps_transport_coeffs_vary_dv:
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g_{}_dv_from_{}_to_{}'.format(U1, U2, U12, T, gamma1,dv_range[0],dv_range[-1])
    title = '$U_1=$'+str(U1)+', $U_2=$'+str(U2)+', $U_{{12}}=$'+str(U12)+', $T=$'+str(T)+', $\\gamma=$'+str(gamma)



    G,S,kappa ,ZT = compute_colormaps_transport_coeffs_vary_dv(v_range,dv_range)
    plot_colormaps_transport_coeffs(v_range,dv_range,G,S,kappa ,ZT,title,'$\\delta v$','colormaps',output)


if Spectral_function_G01:
    v1=-0.3
    v2=-0.5
    EOM = Equations_of_motion(v1, v2, V, T, T , U1, U2, U12, gamma1, gamma2)
    _correlators_ = optimize.root( EOM.Solve_coupled_system, [0,0,0,0,0,0,0,0,0], method='hybr').x
    #n_approxG1 = optimize.root( EOM.dens_approx_G, [EOM.n1,EOM.n2], method='hybr').x
    print('EOM n1 = {}, n2= {}'.format(EOM.n1,EOM.n2))
    print(_correlators_)


if compute_vHxci:
    n1_values = np.zeros((len(v1_range), len(v2_range))) # For z1 to z6
    n2_values = np.zeros((len(v1_range), len(v2_range)))  # For z1 to z6
    vHxc1 = np.zeros((len(v1_range), len(v2_range)))
    vHxc2 = np.zeros((len(v1_range), len(v2_range)))
    for i, v1 in enumerate(v1_range):
        print(v1)
        for j, v2 in enumerate(v2_range):
            EOM = Equations_of_motion(v1, v2, V, T + 0.5*delta_T, T -0.5*delta_T, U1, U2, U12,gamma1, gamma2)
            n1_int = EOM.n1
            n2_int = EOM.n2
            EOM = Equations_of_motion(0., 0., V, T + 0.5*delta_T, T -0.5*delta_T, 0., 0., 0., gamma1, gamma2)
            vs = optimize.root(EOM.solve_vs,[v1,v2], args =(EOM.n1,EOM.n2) ).x
            n1_values[i, j] = 2*n1_int
            n2_values[i, j] = 2*n2_int
            vHxc1[i, j] = v1-vs[0]
            vHxc2[i, j] = v2-vs[1]
    output = 'U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dT_{}_V_{}'.format(U1, U2, U12, T, gamma1, gamma2,delta_T,V)
    plot_colormap_vHxci(n1_values, n2_values, vHxc1,vHxc2,'','vHxc',output)



print('DONE')

