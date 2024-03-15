import numpy as np
from scipy import integrate
from tqdm import tqdm

from inputs import *
from EOM import *
from GCE import *
from generate_plots import *

#######################################################
#######################################################
Uave = 0.5*(U1+U2)

if U12<=U1 and U12<=U2: ##REGIME 1
    print('REGIME 1')
elif U12>U1 and Uave>U12: #REGIME 2
    print('REGIME 2')
elif Uave>=U1 and U12>=Uave: #REGIME 3
    print('REGIME 3')

list_file = []

if calculate_densities_vs_gate:
    imp = GCE_Impurity(U1,U2,U12)
    for veq in tqdm(v_range):
        imp.set_gate( veq+delta_v/2,veq-delta_v/2)
        dens_GCE = imp.comp_dens(beta)
        EOM = Equations_of_motion( veq+delta_v/2,veq-delta_v/2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
        I_EOM0,Q_EOM0 = EOM.Current(EOM.n1,EOM.n2)
        if delta_T!=0.:efficiency_EOM0, output_power_EOM0 = EOM.efficiency(I_EOM0,Q_EOM0)
        else: efficiency_EOM0, output_power_EOM0 = 0.,0.,0.
        if compute_EOMH:
            dens_EOMH = optimize.root( EOM.solve_H, [0.5*dens_GCE[0],0.5*dens_GCE[1]], method='hybr', tol=1.e-12).x
            I_EOMH,Q_EOMH = EOM.Current_EOMH(dens_EOMH[0],dens_EOMH[1])
            if delta_T!=0.:efficiency_EOMH, output_power_EOMH = EOM.efficiency(I_EOMH,Q_EOMH)
            else: efficiency_EOMH, output_power_EOMH,c_eff_EOMH = 0.,0.,0.
            list_file.append((veq,dens_GCE[0],dens_GCE[1],2*dens_EOMH[0],2*dens_EOMH[1],2*EOM.n1,2*EOM.n2,\
                I_EOMH,I_EOM0,Q_EOMH,Q_EOM0,efficiency_EOMH, output_power_EOMH,efficiency_EOM0, output_power_EOM0))
        else: list_file.append((veq,dens_GCE[0],dens_GCE[1],2*EOM.n1,2*EOM.n2, I_EOM0,Q_EOM0,efficiency_EOM0, output_power_EOM0))
    plot_calculate_densities(list_file)
    

###############################################
###############################################

if calculate_spectral_function:
    EOM = Equations_of_motion(v1,v2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
    n1 = EOM.n1
    n2 = EOM.n2
    if compute_EOMH:
        dens_EOMH = optimize.root( EOM.solve_H, [0.5,0.5], method='hybr').x
    for w in w_range: 
        A_1 = EOM.Ai(w,v1, v2, U1, U2, gamma1, gamma2,EOM.n1,EOM.n2)
        A_2 = EOM.Ai(w,v2, v1, U2, U1, gamma2, gamma1,EOM.n1,EOM.n2)
        A0 = (A_1+A_2)
        if compute_EOMH:
            AH_1 = EOM.Ai_Hartree(w, v1, U1, gamma1,dens_EOMH[0],dens_EOMH[1])
            AH_2 = EOM.Ai_Hartree(w, v2, U2, gamma2,dens_EOMH[1],dens_EOMH[0])
            AH = (AH_1+AH_2)
            list_file.append((w,A0,AH))
        else: list_file.append((w,A0))
    norma_A1_O, err = integrate.quad(EOM.Ai, -np.infty, np.infty,  args =(v1,v2, U1,U2, gamma1,gamma2,EOM.n1,EOM.n2))
    if compute_EOMH:
        norma_A1_H, err = integrate.quad(EOM.Ai_Hartree, -np.infty, np.infty,  args =(v1, U1, gamma1,dens_EOMH[0],dens_EOMH[1]))
        print('|A1|_H = {}, |A1|_0 = {}'.format(norma_A1_H,norma_A1_O))
    else: print('|A1|_0 = {}'.format(norma_A1_O))
    plot_calculate_spectral_function(list_file)



###############################################
###############################################


if colormap_spectral_function:
    A = np.zeros((len(U12_range), len(w_range)))
    for i, U12 in enumerate(tqdm(U12_range)):
        v1 = -0.5*U1-U12
        v2 = -0.5*U2-U12
        for j, w in enumerate(w_range):
            EOM = Equations_of_motion(v1,v2,V,T+0.5*delta_T,T-0.5*delta_T,U1,U2,U12,gamma1,gamma2)
            A_1 = EOM.Ai(w,v1, v2, U1, U2, gamma1, gamma2,EOM.n1,EOM.n2)
            A_2 = EOM.Ai(w,v2, v1, U2, U1, gamma2, gamma1,EOM.n2,EOM.n1)
            A[i, j] = (A_1+A_2)
    plot_colormap_spectral_function(w_range, U12_range, A)


###############################################
###############################################
if colormaps_currents:
    dens_0 =  np.zeros((len(V_range), len(v_range)))
    dens_1 =  np.zeros((len(V_range), len(v_range)))
    I_EOM0 =  np.zeros((len(V_range), len(v_range)))
    Q_EOM0 =  np.zeros((len(V_range), len(v_range)))
    for i, V in enumerate(tqdm(V_range)):
        for j, v in enumerate(v_range):
            EOM = Equations_of_motion(v+delta_v/2,v-delta_v/2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            dens = EOM.densities()
            dens_0[i, j] = 2*EOM.n1
            dens_1[i, j] = 2*EOM.n2
            I, Q= EOM.Current(EOM.n1,EOM.n2)
            I_EOM0[i, j] = I
            Q_EOM0[i, j] = Q
    plot_colormaps_currents(v_range,V_range,dens_0,dens_1,I_EOM0,Q_EOM0)


###############################################
###############################################

if colormap_efficiency:
    eficciency_norm_EOM0 = np.zeros((len(V_range), len(v_range)))  # Initialize the efficiency array
    for i, V in enumerate(tqdm(V_range)):
        for j, v in enumerate(v_range):
            v1 = v + delta_v
            v2 = v - delta_v
            EOM = Equations_of_motion(v1, v2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            I, Q = EOM.Current(EOM.n1, EOM.n2)
            eficciency_norm_EOM0[i, j], _ = EOM.efficiency(I, Q)
    plot_colormap_efficiency(v_range,V_range,eficciency_norm_EOM0)

###############################################
###############################################

if colormap_residues_stability_diagram:
    Z = np.zeros((len(v1_range), len(v2_range), 6))  # For z1 to z6
    I_values = np.zeros((len(v1_range), len(v2_range)))
    Q_values = np.zeros((len(v1_range), len(v2_range)))
    nsum = np.zeros((len(v1_range), len(v2_range)))
    for i, v1 in enumerate(tqdm(v1_range)):
        for j, v2 in enumerate(v2_range):
            EOM = Equations_of_motion(v1, v2, V, T + 0.5 * delta_T, T - 0.5 * delta_T, U1, U2, U12, gamma1, gamma2)
            nsum[i,j] = (2*EOM.n1+5*EOM.n2)
            Z[i, j, :] = EOM.numerators(v2, v1, U2, U1, gamma2, gamma1, EOM.n2, EOM.n1)
            I, Q = EOM.Current(EOM.n1, EOM.n2)
            I_values[i, j] = I*100
            Q_values[i, j] = Q*100
    plot_colormap_residues_stability_diagram(v1_range, v2_range, I_values, Q_values,Z, nsum)


###############################################
###############################################

print('DONE')
