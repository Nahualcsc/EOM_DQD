import numpy as np

T=0.1
gamma = 0.2
V= 0.
delta_T = 0.0
U1 = 1. 
U2 =  1.
U12 = 1.



save_data = False
compute_EOMH= False

## uncoment the lines under the calculation to perform 

calculate_densities_vs_gate = False
#v_range = np.linspace(-9,2,100)
#delta_v = 0.

calculate_spectral_function = False
#w_range = np.linspace(-4,4,300)
#v1 = -0.5*U1-U12
#v2 = -0.5*U2-U12

calculate_transport_coeffs = False
#v_range = np.linspace(-4,2,250)
#delta_v = 0.

colormap_spectral_function = False 
#w_range = np.linspace(-3,3,300)
#U12_range = np.linspace(0, 4., 250) 
#delta_v = 0.
## at ph , to change v1,v2 do it at line 80 of main.py


colormaps_transport_coeffs_vary_U12 = False
#delta_v = 0.
#v_range = np.linspace(-5,2,100)
#U12_range = np.linspace(0, 4., 250) 

colormaps_transport_coeffs_vary_dv = False
#v_range = np.linspace(-5,2,100)
#dv_range = np.linspace(-2, 2., 250) 

movie_spectral_function_vary_dv = False
#dv_range = np.linspace(0 , 2, 30)

movie_spectral_function_vary_V = False
#delta_v = 0.
#V_range = np.linspace(0 , 2, 20)

movie_spectral_function_vary_dv = False
#dv_range = np.linspace(0 , 2, 30)


colormaps_currents = False
#delta_v = 0.
#v_range = np.linspace(-9, 2, 100)
#V_range = np.linspace(-4, 4, 100)

movie_currents_vary_U12 = False
#U12_range = np.linspace(0., 3, 40)

colormap_efficiency = False
#delta_v = 0.
#v_range = np.linspace(-4 , 1, 100)
#V_range = np.linspace(-0.2, 0.2, 100)

movie_efficiency_vary_DT = False
#v_range = np.linspace(-5.5 , 1.5, 150)
#V_range = np.linspace(-0.6, 0.6, 150)
#DT_range = np.linspace(0.01, T-0.02, 40)

colormap_SD = False
#v1_range = np.linspace(-6, 1., 100)
#v2_range = np.linspace(-6, 1., 100)

movie_SD_vary_V = False
#V_range = np.linspace(0., 1.5, 40)

movie_SD_vary_U12= False
#U12_range = np.linspace(0., 3, 40)

movie_SD_vary_DT = False
#DT_range = np.linspace(0.01, T-0.02, 40)


Spectral_function_G01 = False



compute_vHxci = False
#v1_range = np.linspace(-6, 1., 100)
#v2_range = np.linspace(-6, 1., 100)



##########
beta= 1./T
gamma1 = gamma 
gamma2 = gamma
