import numpy as np

T=0.1
gamma = 0.1
V= 0.04
delta_T = 0.01
U1 = 2. 
U2 =  3.
U12 = 1.

## uncoment the lines under the calculation to perform 

calculate_densities_vs_gate = False
#v_range = np.linspace(-2,1,30)
#delta_v = 0.

calculate_spectral_function = True
#w_range = np.linspace(-4,4,100)
#v1 = -0.5*U1-U12
#v2 = -0.5*U2-U12

colormap_spectral_function = False 
#w_range = np.linspace(-4,4,200)
#U12_range = np.linspace(0, 8., 200) 
### v1,v2 at ph, can be changed in line 83 main.py

colormaps_currents = False
#delta_v = 0.
#v_range = np.linspace(-5, 2, 100)
#V_range = np.linspace(-4, 4, 100)

colormap_efficiency = False
#delta_v = 0.
#v_range = np.linspace(-4 , 1, 100)
#V_range = np.linspace(-0.2, 0.2, 100)

colormap_residues_stability_diagram = False
#v1_range = np.linspace(-6, 1., 100)
#v2_range = np.linspace(-6, 1., 100)


compute_EOMH= False
#implemented for: calculate_densities_vs_gate, calculate_spectral_function
## T < 0.1 might can be incorrect (EOMH)

##########
beta= 1./T
gamma1 = gamma2 = gamma 
