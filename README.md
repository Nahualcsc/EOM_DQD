# EOM_DQD
Equation of Motion Approach for a parallel double quantum dot system symmetrically attached to two non-interacting electron reservoirs.  The Hamiltonian of the system has the form

![image](https://github.com/Nahualcsc/EOM_DQD/assets/33580847/eda646d0-275a-4a1f-b317-a07837aed386)

where $\alpha=1,2$ refers to the quantum dots and $i=L,R$ to the electron reservoirs. 

The spectral function of the system $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$, the local occupations $n_{\alpha}$, the density correlators $\braket{n_\alpha\ldots n_{\alpha}}$, as well as the  charge  and heat currents I and Q are analytically calculated for a given out-of-equilibrium working configuration set by the parameters $T,\gamma$, the external gates $\varepsilon_\alpha$, the bias voltage $V$, the thermal gradient $\Delta T=T_L-R_R$ and the Coulombic interactions $U_\alpha$ and $U_{12}$.

## Parameters, Interactions, and calculation setup
There are currently six possible calculations implemented
1. ### calculate_densities_vs_gate:
   It computes the local occupations $n_\alpha$ as a function of the external gate $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$. The steady-state charge $I$ and heat $Q$ currents are also computed in a non-equilibrium setup. For finite thermal gradient $\delta T$, the thermal efficiency (normalized over the Carnot efficiency) is also computed.
2. ### calculate_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the frequency $\omega$.
3. ### colormap_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the inter-Coulomb interaction $U_{12}$ and  frequency $\omega$.
6. colormaps_currents
7. colormap_efficiency_vary_U
8. colormap_correlators_stability_diagram
## Program usage
The program can be run by typing in a linux shell: ```python3 inputs.py ```

## System requirements
Python3 is required and the following Python packages:
- numpy
- scipy
- maplotlib
- tqdm
  
