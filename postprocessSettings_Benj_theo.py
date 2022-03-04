###################################################################
# SIMULATION INFOS
###################################################################

# # L=16*Pi*H / epsilon=0.2
# current_path = 'BL_Vortices_Eps0d2_Thet0_L'
# mean_path = 'Preliminary_Laminar_BLFlow_STEP4'
# Re = 1000
# L1 = 16*np.pi
# L3 = 8*np.pi

# number_iteration = 170
# preliminary_iterations = 85

# # L=32*Pi*H / epsilon=0.2 / theta=45°
# current_path = 'BL_Vortices_Eps0d2_Thet45_2L'
# mean_path = 'BL_Re1000_2L'
# Re = 1000
# L1 = 32*np.pi
# L3 = 8*np.pi

# number_iteration = 309np.nan_to_num(
# # number_iteration = 254 # 164 + 90
# preliminary_iterations = 163

# # L=32*Pi*H / epsilon=0.4 / theta=0°
# current_path = 'BL_Vortices_Eps0d4_Thet0_2L'
# mean_path = 'BL_Re1000_2L'
# Re = 1000
# L1 = 32*np.pi
# L3 = 8*np.pi

# number_iteration = 308
# # number_iteration = 254 # 164 + 90
# preliminary_iterations = 163

# # L=32*Pi*H / epsilon=0.4 / theta=45°
# current_path = 'BL_Vortices_Eps0d4_Thet45_2L'
# mean_path = 'BL_Re1000_2L'
# Re = 1000
# L1 = 32*np.pi
# L3 = 8*np.pi

# number_iteration = 306
# # number_iteration = 254 # 164 + 90
# preliminary_iterations = 163

# # L=32*Pi*H / epsilon=0.00001 / theta=0°
# current_path = 'PL_Vortices_SmallEps'
# mean_path = 'Preliminary_Laminar_PFlow'
# Re = 3000
# Pr = 1
# deltaT = 0.5
# L1 = 1.5*32*np.pi
# L3 = 8*np.pi
# symmetry = True

# number_iteration = 160
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.07 / theta=0°
# current_path = 'PL_Vortices_ModerateEps'
# mean_path = 'Preliminary_Laminar_PFlow'
# Re = 3000
# Pr = 1
# deltaT = 0.5
# L1 = 1.5*32*np.pi
# L3 = 8*np.pi
# symmetry = True

# number_iteration = 160
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Dir'
# Re = 3000
# Pr = 1
# deltaT = 0.5
# L1 = 1.5*32*np.pi
# L3 = 8*np.pi
# symmetry = True

# number_iteration = 160
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps_Ka_coarse'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Ka'
# Re = 3000
# Pr = 1
# deltaT = 0
# L1 = 32*np.pi
# L3 = 8*np.pi
# symmetry = True
# delta_iteration = 5

# number_iteration = 175
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps_Ka_NW'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Ka'
# Re = 3000
# Pr = 1
# deltaT = 0
# L1 = 32*np.pi
# L3 = 8*np.pi
# symmetry = False
# delta_iteration = 5

# number_iteration = 175
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps_Prel'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Dir'
# Re = 3000
# Pr = 1
# deltaT = 0.5
# L1 = 32*np.pi
# L3 = 8*np.pi
# symmetry = True

# number_iteration = 50
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps_fine'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Dir'
# Re = 3000
# Pr = 1
# deltaT = 0.5
# L1 = 32*np.pi
# L3 = 8*np.pi
# symmetry = True
# delta_iteration = 5

# number_iteration = 175
# preliminary_iterations = 10

# # L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps_Ka'
# mean_path = 'Preliminary_Laminar_PFlow'
# theo_path = 'TF_ReTau_150_PS_Dir'
# Re = 3000
# Pr = 1
# deltaT = 0.1
# L1 = 32*np.pi
# L3 = 8*np.pi
# symmetry = True
# delta_iteration = 5

# number_iteration = 175
# preliminary_iterations = 10


# L=4*Pi*H / epsilon=0.14 / theta=0°
current_path = 'ReTau_150_PS_Dir'
mean_path = 'ReTau_150_PS_Dir'
theo_path = ''
Re = 3000
Pr = 1
deltaT = 0.5
L1 = 4*np.pi
L3 = (4/3)*np.pi
symmetry = True
delta_iteration = 5

number_iteration = 170
preliminary_iterations = 90

root = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/'
root_postproc = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/PostProcessing/'

###################################################################
# CLASSIC PLOT INFOS
###################################################################

custom_t = False
# rangeT = [-6,-1,4,9,14,19,24,29,34,39,44,49,54,59]
# rangeT = [94,99,104,109,114,119,124,129,134,139,144,149,154,159,164]
# rangeT = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]
# rangeT = [-6,-1]
# rangeT = [139, 144, 149, 154, 159, 164]
rangeT = [165]

# 3D plotting
plot_U_fluctuations = False
plot_V_fluctuations = False
plot_W_fluctuations = False

plot_U_instantaneous = False
plot_V_instantaneous = False
plot_W_instantaneous = False

plot_T_fluctuations = False
plot_T_instantaneous = False

plot_x_equal_to = False
plot_y_equal_to = 1
plot_z_equal_to = False #3.5*np.pi

shift = 77.5

plot_spot = False
plot_thermal_spot = False
generate_video_spot = False

plot_lambda2 = False
generate_video_lambda2 = False

dwdx_plot = False
generate_dwdx_video = False

colored_contours = False
contours = True

compute_lambda2 = False
slice_lambda2 = 64

# 2D plotting
plot_wall_shear_stress_contours = False

plot_wall_temperature_flux_contours = False

# 1D stats generation
compute_velocity_stats = True
compute_transport_equation_terms = True
compute_reynolds_stresses = True
compute_vorticity = True
compute_temperature_stats = True
compute_heat_flux = True
compute_temperature_velocity_correlations = True
compute_temperature_budgets = True

compute_E_tot_evolution = False

# 1D stats plotting
plot_velocity_stats = True
plot_transport_equation_terms = True
plot_reynolds_stresses = True
plot_vorticity = True
plot_temperature_stats = True
plot_temperature_velocity_correlations = True
plot_temperature_budgets = True

plot_E_tot_evolution = False
plot_Re_tau_evolution = False
plot_T_tau_evolution = False
plot_Nu_evolution = False


# choose spot
total_spot = False
heart_spot = False
wave_packet_spot = False
# for the domain
total_domain = True

x_bounds = [24*np.pi, 32*np.pi]
y_bounds = [0, 2]
z_bounds = [0, 8*np.pi]

###################################################################
# Henningson PLOT INFOS
###################################################################

# plot temporal evolution like Henningson Fig 2
plot_TEHL_U_fluctuations = False
plot_TEHL_V_fluctuations = False
plot_TEHL_W_fluctuations = False
plot_TEHL_T_fluctuations = False

plot_TEHL_wall_shear_stress_contours = False

plot_TEHL_wall_temperature_flux_contours = False

# 1D stats plotting
plot_TEHL_velocity_stats = False
plot_TEHL_transport_equation_terms = False
plot_TEHL_reynolds_stresses = False
plot_TEHL_temperature_stats = False
plot_TEHL_temperature_velocity_correlations = False

# # L=16*Pi*H / epsilon=0.2
# # current_path = 'Laminar_BLFlow_CounterRotatingVortices'
# # Re = 1000
# time_values = [10,20,30,40]
# y_limits = [[5,20], [5,20], [5,20], [5,20]]
# x_limits = [[0,30], [10,40], [20,50], [30,60]]

# # L=32*Pi*H / epsilon=0.2 / theta=45°
# # current_path = 'BL_Vortices_Eps0d2_Thet45_2L'
# # Re = 1000
# time_values = [10,20,30,40]
# y_limits = [[5,20], [5,20], [5,20], [5,20]]
# x_limits = [[0,30], [10,40], [20,50], [30,60]]

# # L=32*Pi*H / epsilon=0.4 / theta=0°
# # current_path = 'BL_Vortices_Eps0d4_Thet0_2L'
# # Re = 1000
# time_values = [10,20,30,40]
# y_limits = [[5,20], [5,20], [5,20], [5,20]]
# x_limits = [[0,30], [10,40], [20,50], [30,60]]

# # L=32*Pi*H / epsilon=0.4 / theta=45°
# # current_path = 'BL_Vortices_Eps0d4_Thet45_2L'
# # Re = 1000
# time_values = [10,20,40,80]
# y_limits = [[0,25], [0,25], [0,25], [0,25]]
# x_limits = [[0,30], [20,60], [50,110], [90,115]]

# # L=32*Pi*H / epsilon=0.00001 / theta=0°
# # current_path = 'PL_Vortices_SmallEps'
# # Re = 3000
# time_values = [19,39,79,119]
# y_limits = [[0,25], [0,25], [0,25], [0,25]]
# x_limits = [[0,30], [10,40], [20,70], [30,100]]

# L=32*Pi*H / epsilon=0.07 / theta=0°
# current_path = 'PL_Vortices_SmallEps'
# Re = 3000
# time_values = [9,19,29,39]
# y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]
# x_limits = [[0,20], [5,30], [10,40], [10,50]]
# # wall
# time_values = [19,39,79,119]
# y_limits = [[0,25], [0,25], [0,25], [0,25]]
# x_limits = [[0,30], [10,40], [20,70], [30,100]]
# current_path = 'PL_Vortices_ModerateEps'
# Re = 3000
# time_values = [19,29,49,89]
# y_limits = [[0,25], [0,25], [0,25], [0,25]]
# x_limits = [[0,30], [5,40], [10,60], [20,110]]

# # L=32*Pi*H / epsilon=0.07 / theta=0°
# # current_path = 'PL_Vortices_ModerateEps'
# # Re = 3000
# time_values = [9,19,29,39]
# y_limits = [[-10,10], [-10,10], [-10,10], [-10,10]]
# x_limits = [[0,20], [5,30], [10,35], [15,45]]
# time_values = [19,29,49,89]
# y_limits = [[0,25], [0,25], [0,25], [0,25]]
# x_limits = [[0,30], [5,40], [10,60], [20,110]]
# # wall
# time_values = [-6,14,29,44]
# y_limits = [[-10,10], [-10,10], [-10,10], [-10,10]]
# x_limits = [[0,30], [10,40], [20,70], [30,100]]

# time_values = [19,39,79,119]
# y_limits = [[-10,10], [-10,10], [-10,10], [-10,10]]
# y_limits = [[0,2], [0,2], [0,2], [0,2]]
# x_limits = [[0,30], [20,60], [30,90], [50,120]]

# L=32*Pi*H / epsilon=0.14 / theta=0°
# current_path = 'PL_Vortices_LargeEps'
# Re = 3000
# time_values = [50,70,90,110]
# y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]
# x_limits = [[10,40], [15,45], [20,55], [25,65]]
# time_values = [19,29,49,89]
# y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]
# x_limits = [[0,30], [5,30], [10,45], [25,75]]
# DEBUG
# time_values = [99,119,139,149]
# y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]
# x_limits = [[40,100], [40,120], [50,120], [60,120]]
# time_values = [39,49,59,69]
# y_limits = [[-10,10], [-10,10], [-10,10], [-10,10]]
# x_limits = [[10,45], [15,50], [20,60], [25,70]]
time_values = [89,109,129,149]
y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]
x_limits = [[10,45], [15,50], [20,60], [25,70]]
