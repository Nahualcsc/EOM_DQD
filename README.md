# EOM_DQD
Equation of Motion Approach for a parallel double quantum dot system symmetrically attached to two non-interacting electron reservoirs.  The Hamiltonian of the system has the form

![image](https://github.com/Nahualcsc/EOM_DQD/assets/33580847/eda646d0-275a-4a1f-b317-a07837aed386)

where $\alpha=1,2$ refers to the quantum dots and $i=L,R$ to the electron reservoirs. 

The spectral function of the system $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$, the local occupations $n_{\alpha}$, the density correlators $\braket{n_\alpha\ldots n_{\alpha}}$, as well as the  charge  and heat currents I and Q are analytically calculated for a given out-of-equilibrium working configuration set by the parameters $T,\gamma$, the external gates $\varepsilon_\alpha$, the bias voltage $V$, the thermal gradient $\Delta T=T_L-T_R$ and the Coulombic interactions $U_\alpha$ and $U_{12}$.

## Parameters, Interactions, and calculation setup
The parameters and interactions to be always set up in ```inputs.py``` are: ```T, delta_T, gamma, V, U1, U2, U12```. 
There are currently six possible calculations implemented
1. ### calculate_densities_vs_gate:
   It computes the local occupations $n_\alpha$ as a function of the external gate $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$. The steady-state charge $I$ and heat $Q$ currents are also computed in a non-equilibrium setup. For finite thermal gradient $\delta T$, the thermal efficiency $\eta/\eta_C$(normalized over the Carnot efficiency) is also computed.
   Extra parameters to define: ```v_range, delta_v```.
2. ### calculate_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the frequency $\omega$.
   Extra parameters to define: ```w_range, v1, v2```.
3. ### colormap_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the inter-Coulomb interaction $U_{12}$ and  frequency $\omega$.
   Extra parameters to define: ```w_range, U12_range```.
   Implemented at the particle-hole symmetry point $\varepsilon_\alpha=-U_\alpha/2-U_{12}$, to set other gate levels go to line 78 of main.py
4. ### colormaps_currents:
   The local occupations, as well as the steady-state charge $I$ and heat $Q$ currents, are computed as a function of the external gate $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$ and the bias voltage $V$. 
   Extra parameters to define: ```delta_v, v_range, V_range```.
5. ### colormap_efficiency_vary_U:
   The thermal efficiency  $\eta/\eta_C$ is computed as a function of the $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$ and the bias voltage $V$.
   Extra parameters to define: ```delta_v, v_range, V_range```.
6. ### colormap_residues_stability_diagram:
   The residues (numerators of the spectral function), the stability diagram ($n_1+5*n_2$), and the steady-state charge $I$ and heat $Q$ currents as a function of the external gates $\varepsilon_1$ and  $\varepsilon_2$.
    Extra parameters to define: ```v1_range, v2_range```.

For ```calculate_densities_vs_gate``` and ```calculate_spectral_function``` the calculation can be compared against the Hartree EOM (EOMH) result. To activate it: ```compute_EOMH= False```. Note that EOMH takes more time to evaluate.

## Program usage
The program can be run by typing in a linux shell: ```python3 inputs.py ```

## System requirements
Python3 is required and the following Python packages:
- numpy
- scipy
- maplotlib
- tqdm
  
