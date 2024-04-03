# EOM_DQD
Equation of Motion Approach for a parallel double quantum dot system symmetrically attached to two non-interacting electronic reservoirs.  The Hamiltonian of the system has the form

![image](https://github.com/Nahualcsc/EOM_DQD/assets/33580847/eda646d0-275a-4a1f-b317-a07837aed386)

where $\alpha=1,2$ refers to the quantum dots and $i=L,R$ to the electron reservoirs. 

The spectral function of the system $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$, the local occupations $n_{\alpha}$, the density correlators $\braket{n_\alpha\ldots n_{\alpha}}$, as well as the  charge  and heat currents I and Q are analytically calculated for a given out-of-equilibrium working configuration set by the parameters $T,\gamma$, the external gates $\varepsilon_\alpha$, the bias voltage $V$, the thermal gradient $\Delta T=T_L-T_R$ and the Coulombic interactions $U_\alpha$ and $U_{12}$.

## Parameters, Interactions, and calculation setup
The parameters and interactions to be always set up in ```inputs.py``` are: ```T, delta_T, gamma, V, U1, U2, U12```. 
There are currently six possible calculations implemented and another seven animations related. In the following $\varepsilon_\alpha=v_\alpha$.
1. ### calculate_densities_vs_gate:
   The local occupations $n_\alpha$ are computed as a function of the external gate $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$. The steady-state charge $I$ and heat $Q$ currents are also computed in a non-equilibrium setup. For finite thermal gradient $\Delta T$, the thermal efficiency $\eta/\eta_C$ (normalized over the Carnot efficiency) is also computed.
   
   Extra parameters to define: ```v_range, delta_v```.
2. ### calculate_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the frequency $\omega$.
   
   Extra parameters to define: ```w_range, v1, v2```.
3. ### colormap_spectral_function:
   The total spectral function $A=-\sum_{\alpha}\text{Im}\[G_\alpha(\omega)\]/\pi$ is calculated as a function of the inter-Coulomb interaction $U_{12}$ and  frequency $\omega$.
   
   Extra parameters to define: ```w_range, U12_range```.
   
   Implemented at the particle-hole symmetry point $\varepsilon_\alpha=-U_\alpha/2-U_{12}$, to set other gate levels go to line 80 of main.py
4. ### colormaps_currents:
   The local occupations, as well as the steady-state charge $I$ and heat $Q$ currents, are computed as a function of the external gate $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$ and the bias voltage $V$.
   
   Extra parameters to define: ```delta_v, v_range, V_range```.
5. ### colormap_efficiency:
   The thermal efficiency  $\eta/\eta_C$ is computed as a function of the $\varepsilon=\frac{\varepsilon_1+\varepsilon_2}{2}$ and the bias voltage $V$.
   
   Extra parameters to define: ```delta_v, v_range, V_range```.
   
6. ### colormap_SD:
   The residues (numerators of the Green's function), the stability diagram ($n_1+5n_2$), and the steady-state charge $I$ and heat $Q$ currents as a function of the external gates $\varepsilon_1$ and  $\varepsilon_2$.
   
    Extra parameters to define: ```v1_range, v2_range```.

For ```calculate_densities_vs_gate``` and ```calculate_spectral_function``` the calculation can be compared against the Hartree EOM (EOMH) result. To activate it: ```compute_EOMH= True```. Note that EOMH takes more time to evaluate.

 7. ### movie_spectral_function_vary_dv:
    Animation of ```calculate_spectral_function```  for ```dv_range```.
 8. ### movie_spectral_function_vary_V:
     Animation of ```calculate_spectral_function```  for ```dV_range```.
 9. ### movie_currents_vary_U12:
     Animation of ```colormaps_currents``` for ```U12_range```.
 9. ### movie_SD_vary_V:
     Animation of ```colormap_SD``` for ```V_range```.
 11. ### movie_SD_vary_U12:
     Animation of ```colormap_SD``` for ```U12_range```.
 12. ### movie_SD_vary_DT:
      Animation of ```colormap_SD``` for ```DT_range```.
 13. ### movie_efficiency_vary_DT:
     Animation of ```colormap_efficiency``` for ```DT_range```.

## Examples
 1. $U_1=1,U_2=2,T=0.1,\gamma_1=\gamma_2=10^{-3},V=0,\Delta T=10^{-2}$ and $U_{12}\in (0,3)$:
  <img src="https://github.com/Nahualcsc/EOM_DQD/assets/33580847/99e7f1e9-9618-43d9-85d3-764b4c62827a" width="250" alt="(I, Q) from movie_SD_vary_U12">

2. $U_1=1,U_2=2,T=0.1,\gamma_1=\gamma_2=10^{-3},V=0,\Delta T=10^{-2}$ and $U_{12}\in (0,3)$:
  <img src="https://github.com/Nahualcsc/EOM_DQD/assets/33580847/24d2124a-f7e9-4d97-bb88-161c4884e55d" width="250" alt="Stability diagram from movie_SD_vary_U12">

3. $U_1=2,U_2=3,T=0.1,\gamma_1=\gamma_2=10^{-1},\delta v=0,\Delta T=0$ and $U_{12}\in (0,3)$:
  <img src="https://github.com/Nahualcsc/EOM_DQD/assets/33580847/0482499e-1861-4852-bd57-134b0202e457" width="250" alt=" (N,I,Q) from movie_currents_vary_U12">
  
4. $U_1=2,U_2=3,T=0.1,\gamma_1=\gamma_2=10^{-1},\delta v=0,\Delta T=0$ and $U_{12}\in (0,3)$:
  <img src="https://github.com/Nahualcsc/EOM_DQD/assets/33580847/5f82cae8-daed-4dde-a3d7-873bccf824dd" width="250" alt=" (n1,n2) from movie_currents_vary_U12">

5. $U_1=2,U_2=3,T=0.1,\gamma_1=\gamma_2=0,\delta v=0,\Delta T=0, v_\alpha = -U_\alpha/2-U_{12}\pm\delta_v/2$ and $\delta v\in (0,2)$:
  <img src="https://github.com/Nahualcsc/EOM_DQD/assets/33580847/4f0468f1-fda3-420a-a378-3e318067fb63" width="250" alt="Spectral function from movie_spectral_function_vary_dv">



## Program usage
The program can be run by typing in a linux shell: ```python3 main.py ```.

The ```inputs.py``` file contain the neccessary commented flags and inputs to set it up. 

## System requirements
Python3 is required and the following Python packages:
- numpy
- scipy
- maplotlib
- tqdm
- os








