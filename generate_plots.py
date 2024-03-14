import numpy as np
from scipy import integrate
from tqdm import tqdm
from inputs import *
from EOM import *
from GCE import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os

#######################################################
#######################################################

def  plot_calculate_densities(list_file):
    if not os.path.exists('ni'): os.makedirs('ni')
    if not os.path.exists('Currents'): os.makedirs('Currents')
    if not os.path.exists('Efficiency'): os.makedirs('Efficiency')
    fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharex=True)
    if compute_EOMH:    
        veq,n1GCE,n2GCE,n1EOMH,n2EOMH,n1EOM0,n2EOM0,I_EOMH,I_EOM0,Q_EOMH,Q_EOM0,\
        efficiency_EOMH, output_power_EOMH,efficiency_EOM0, output_power_EOM0 = zip(*list_file)
        axs[0].plot(veq,n1EOMH, '--', color='red',label = 'EOM H',linewidth=1.5)
        axs[1].plot(veq,n2EOMH, '-', color='red',label = 'EOM H',linewidth=1.5)
    else:    
        veq,n1GCE,n2GCE,n1EOM0,n2EOM0,I_EOM0,Q_EOM0,efficiency_EOM0, output_power_EOM0 = zip(*list_file)
    axs[0].plot(veq,n1GCE, '-', color='green',label = 'GCE',linewidth=1.5)
    axs[0].plot(veq,n1EOM0, '--', color='blue',label = 'EOM 0',linewidth=1.5)
    axs[1].plot(veq,n2GCE, '-', color='green',label = 'GCE',linewidth=1.5)
    axs[1].plot(veq,n2EOM0, '--', color='blue',label = 'EOM 0',linewidth=1.5)
    axs[0].set_ylabel('$n_1$', fontsize=16)
    axs[1].set_ylabel('$n_2$', fontsize=16)
    axs[1].set_xlabel('$v$', fontsize=16)
    fig.suptitle('$U1$ = {}, $U2$ = {}, $U12$ = {}, $T$ = {}, $\\gamma_1$ = {}, $\\gamma_2$ = {}, $\\delta v$ = {}'.format(U1,U2,U12,T,gamma1,gamma2,delta_v))
    fig.tight_layout()
    axs[0].legend(loc='upper left')
    fig.savefig('./ni/n_V_{}_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_dv_{}.png'.format(V,U1,U2,U12,T,gamma1,gamma2,delta_v), dpi=300)

    if V!=0. or delta_T!=0.:
        fig2, axs2 = plt.subplots(2, 1, figsize=(8, 12), sharex=True)
        if compute_EOMH:    
            axs2[0].plot(veq,I_EOMH, '-', color='red',label = 'EOM H',linewidth=1.5)
            axs2[1].plot(veq,Q_EOMH, '-', color='red',label = 'EOM H',linewidth=1.5)
        axs2[0].plot(veq,I_EOM0, '--', color='blue',label = 'EOM 0',linewidth=1.5)       
        axs2[1].plot(veq,Q_EOM0, '--', color='blue',label = 'EOM 0',linewidth=1.5)    
        axs2[0].set_ylabel('$I/\\gamma$', fontsize=16)
        axs2[1].set_ylabel('$Q/\\gamma$', fontsize=16)
        axs2[1].set_xlabel('$v$', fontsize=16)
        fig2.suptitle('$V$ = {}, $U1$ = {}, $U2$ = {}, $U12$ = {}, $T$ = {}, $\\Delta T$ = {} $\\gamma_1$ = {}, $\\gamma_2$ = {}, $\\delta v$ = {}'.format(V,U1,U2,U12,T,delta_T,gamma1,gamma2,delta_v))
        fig2.tight_layout()
        axs2[0].legend(loc='upper left')
        fig2.savefig('./Currents/IQ_V_{}_U1_{}_U2_{}_U12_{}_T_{}_delta_T_{}_g1_{}_g2_{}_dv_{}.png'.format(V,U1,U2,U12,T,delta_T,gamma1,gamma2,delta_v), dpi=300)
    if delta_T!=0.:
        fig2, axs2 = plt.subplots(2, 1, figsize=(8, 12), sharex=True)
        if compute_EOMH:   
            axs2[0].plot(veq,efficiency_EOMH, '--', color='red',label = 'EOM H',linewidth=1.5)
            axs2[1].plot(veq,output_power_EOMH, '--', color='red',label = 'EOM H',linewidth=1.5)
        axs2[0].plot(veq,efficiency_EOM0, '-', color='blue',label = 'EOM 0',linewidth=1.5)       
        axs2[1].plot(veq,output_power_EOM0, '-', color='blue',label = 'EOM 0',linewidth=1.5)    
        axs2[0].set_ylabel('$\\eta/\\eta_C$', fontsize=16)
        axs2[1].set_ylabel('$P$', fontsize=16)
        axs2[1].set_xlabel('$v$', fontsize=16)
        fig2.suptitle('$V$ = {}, $U1$ = {}, $U2$ = {}, $U12$ = {}, $T$ = {}, $\\Delta T$ = {} $\\gamma_1$ = {}, $\\gamma_2$ = {}, $\\delta v$ = {}'.format(V,U1,U2,U12,T,delta_T,gamma1,gamma2,delta_v))
        fig2.tight_layout()
        axs2[0].legend(loc='upper left')
        fig2.savefig('./Efficiency/eta_V_{}_U1_{}_U2_{}_U12_{}_T_{}_delta_T_{}_g1_{}_g2_{}_dv_{}.png'.format(V,U1,U2,U12,T,delta_T,gamma1,gamma2,delta_v), dpi=300)

###############################################
###############################################


def plot_calculate_spectral_function(list_file):
    if not os.path.exists('Aw'): os.makedirs('Aw')
    fig, axs = plt.subplots(1, 1, figsize=(6, 5), sharex=True)
    if compute_EOMH:
        w,A0,AH= zip(*list_file)
        axs.plot(w,AH, '-', color='blue',label = 'EOM H',linewidth=1.5)
    else: w,A0= zip(*list_file)
    axs.plot(w,A0, '-', color='red',label = 'EOM 0',linewidth=1.5)
    axs.set_ylabel('$A(\\omega)$', fontsize=16)
    axs.set_xlabel('$\\omega$', fontsize=16)
    fig.tight_layout()
    axs.legend(loc='upper left')
    fig.savefig('./Aw/A_v1_{}_v2_{}_V_{}_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_delta_T_{}.png'.format(v1,v2,V,U1,U2,U12,T,gamma1,gamma2,delta_T), dpi=300)



###############################################
###############################################



def plot_colormap_spectral_function(w_range, U12_range, A):
    if not os.path.exists('Aw'): os.makedirs('Aw')
    fig1, ax1 = plt.subplots(1, 1, figsize=(6,5)) 
    c = ax1.pcolormesh(w_range, U12_range, A, shading='auto', cmap='inferno')
    ax1.set_xlabel('$\\omega$')
    ax1.set_ylabel('$U_{12}$')
    fig1.colorbar(c, ax=ax1)
    fig1.tight_layout()
    fig1.savefig('./Aw/A_U1_{}_U2_{}_T_{}_g1_{}_g2_{}_V_{}_dT_{}.png'.format(U1, U2, T, gamma1, gamma2,V,delta_T), dpi=300)



###############################################
###############################################

def plot_colormaps_currents(v_range,V_range,dens_0,dens_1,I_EOM0,Q_EOM0):
    if not os.path.exists('colormaps'): os.makedirs('colormaps')
    fontsize = 20 
    fig, axs = plt.subplots(2, 1, figsize=(8, 16)) 
    c0 = axs[0].imshow(dens_0, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    c1 = axs[1].imshow(dens_1, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    fig.colorbar(c0, ax=axs[0], label='$n_1$').ax.tick_params(labelsize=fontsize)
    fig.colorbar(c1, ax=axs[1], label='$n_2$').ax.tick_params(labelsize=fontsize)
    axs[0].set_title('EOM 0', fontsize=16)
    axs[0].set_xlabel('$\\varepsilon$', fontsize=fontsize)
    axs[0].set_ylabel('V', fontsize=fontsize)
    axs[1].set_xlabel('$\\varepsilon$', fontsize=fontsize)
    axs[1].set_ylabel('V', fontsize=fontsize)
    axs[0].tick_params(axis='y', which='major', labelsize=fontsize)
    axs[1].tick_params(axis='both', which='major', labelsize=fontsize)
    axs[0].text(0.5, 0.05, '$n_1$', transform=axs[0].transAxes, ha='center', va='center', color='white', fontsize=fontsize)
    axs[1].text(0.5, 0.05, '$n_2$', transform=axs[1].transAxes, ha='center', va='center', color='white', fontsize=fontsize)
    fig.tight_layout(pad=1.0)
    fig.savefig('./colormaps/ni_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}.png'.format(U1, U2, U12, T, gamma1, gamma2), dpi=300)
    
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    c1 = axs[0].imshow(dens_0, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    c2 = axs[1].imshow(I_EOM0, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    c3 = axs[2].imshow(Q_EOM0, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    fig.colorbar(c1, ax=axs[0]).ax.tick_params(labelsize=fontsize)
    fig.colorbar(c2, ax=axs[1]).ax.tick_params(labelsize=fontsize)
    fig.colorbar(c3, ax=axs[2]).ax.tick_params(labelsize=fontsize)
    for ax in axs:
        ax.set_xlabel('$\\epsilon$', fontsize=fontsize)
    axs[0].tick_params(axis='both', which='major', labelsize=fontsize)
    axs[1].tick_params(axis='x', which='major', labelsize=fontsize)
    axs[2].tick_params(axis='x', which='major', labelsize=fontsize)
    axs[0].set_ylabel('$V$', fontsize=fontsize)
    axs[0].text(0.5, 0.05, '$N$', transform=axs[0].transAxes, ha='center', va='center', color='white', fontsize=fontsize)
    axs[1].text(0.5, 0.05, '$I/\\gamma$', transform=axs[1].transAxes, ha='center', va='center', color='white', fontsize=fontsize)
    axs[2].text(0.5, 0.05, '$Q/\\gamma$', transform=axs[2].transAxes, ha='center', va='center', color='white', fontsize=fontsize)
    fig.subplots_adjust(wspace=0.15, hspace=0.3) 
    fig.tight_layout(pad=1.0)  
    fig.savefig('./colormaps/n1IQ_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}.png'.format(U1, U2, U12, T, gamma1, gamma2), dpi=300)

###############################################
###############################################

def plot_colormap_efficiency_vary_U(v_range,V_range,eficciency_norm_EOM0):
    if not os.path.exists('colormaps'): os.makedirs('colormaps')
    fig, axs = plt.subplots(1, 1, figsize=(6, 5), sharex=True) 
    c = axs.imshow(eficciency_norm_EOM0, extent=[v_range.min(), v_range.max(), V_range.min(), V_range.max()], origin='lower', aspect='auto', cmap='inferno')
    cb = fig.colorbar(c, ax=axs)
    cb.set_label(label='$\\eta/\\eta_C$', size=20)  
    cb.ax.tick_params(labelsize=18)  
    axs.set_xlabel('$\\epsilon$', fontsize=20)  
    axs.set_ylabel('V', fontsize=20)  
    axs.tick_params(axis='both', which='major', labelsize=18) 
    plt.tight_layout()
    plt.savefig('./colormaps/efficiency_U1_{}_U2_{}_T_{}_g1_{}_g2_{}.png'.format(U1, U2, T, gamma1, gamma2), dpi=300)

###############################################
###############################################

def plot_colormap_correlators_stability_diagram(v1_range, v2_range, I_values, Q_values,Z, nsum):
    if not os.path.exists('colormaps'): os.makedirs('colormaps')
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))  # Adjusted figure size
    plt.subplots_adjust(left=0.1, right=0.92, bottom=0.1, top=0.95, wspace=0.05, hspace=0.05) 
    titles = ['$z_1$', '$z_2$', '$z_3$', '$z_4$', '$z_5$', '$z_6$']
    for i, (ax, z, title) in enumerate(zip(axs.flat, Z.transpose(2, 0, 1), titles)):
        c = ax.pcolormesh(v1_range, v2_range, z, shading='auto', cmap='inferno', vmin=0, vmax=1)
        ax.text(0.5, 0.9, title, transform=ax.transAxes, ha='center', va='center', color='white', fontsize=20)
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(-5,1)
        ax.tick_params(axis='both', which='major', labelsize=20)
        if i >= 3:
            ax.set_xlabel('$\\varepsilon_1$', fontsize=20)
        else:
            ax.set_xticklabels([])
        if i % 3 == 0:
            ax.set_ylabel('$\\varepsilon_2$', fontsize=20)
        else:
            ax.set_yticklabels([])
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])  # Position for the color bar
    fig.colorbar(c, cax=cbar_ax)
    plt.savefig('./colormaps/residues_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_V_{}_dT_{}.png'.format(U1, U2, U12, T, gamma1, gamma2,V,delta_T), dpi=300)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3 ))  # 1 row, 2 cols for I and Q
    plt.subplots_adjust(wspace=0.4)
    for ax, val, title, vmin, vmax in zip(axs, [I_values, Q_values], ['$I/\\gamma\\cdot 10^2$', '$Q/\\gamma\\cdot 10^2$'], [0,-40],[50,2]):
        c = ax.pcolormesh(v1_range, v2_range, val, shading='auto', cmap='inferno')#, vmin=vmin, vmax=vmax)
        ax.text(0.5, 0.9, title, transform=ax.transAxes, ha='center', va='center', color='white', fontsize=15)
        ax.set_xlabel('$\\varepsilon_1$')
        ax.set_ylim(-5, 1)
        if title == '$I/\\gamma\\cdot 10^2$':
            ax.set_ylabel('$\\varepsilon_2$')
        else:
            # For the right figure, remove y-ticks and y-label
            ax.set_yticks([])
            ax.set_ylabel('')
        fig.colorbar(c, ax=ax)
    fig.subplots_adjust(wspace=0.1)  
    fig.tight_layout()
    plt.savefig('./colormaps/IQ_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_V_{}_dT_{}.png'.format(U1, U2, U12, T, gamma1, gamma2,V,delta_T), dpi=300)

    fig1, ax1 = plt.subplots(1, 1, figsize=(6,5)) 
    c = ax1.pcolormesh(v1_range, v2_range, nsum, shading='auto', cmap='inferno', vmin=0, vmax=7)
    ax1.set_xlabel('$v_1$')
    ax1.set_ylabel('$v_2$')
    ax1.set_ylim(-5, 1)
    fig1.colorbar(c, ax=ax)
    fig1.tight_layout()
    fig1.savefig('./colormaps/SD_U1_{}_U2_{}_U12_{}_T_{}_g1_{}_g2_{}_V_{}_dT_{}.png'.format(U1, U2, U12, T, gamma1, gamma2,V,delta_T), dpi=300)


###############################################
###############################################

