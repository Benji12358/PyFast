###################################################################
################    GENERAL SIMULATION INFOS        ###############
###################################################################

# CHANNEL RE_TAU = 125, w*=12
current_path = 'AC_ReTau_125_w12_PSL'                                                         # Name of the simulation to post-process
mean_path = 'AC_ReTau_125_w12_PSL'                                                            # Name of the simulation containing mean fields (by default, put same name as simulation to post-process)
theo_path = ''                                                                               # Name of the simulation containing theoretical results (by default, leave blank)
Re = 1500                                                                                    # Reynolds number of the simulation
Pr = 1                                                                                       # Prandtl number of the simulation
deltaT = 0.5                                                                                 # Temperature step between the lower and the upper wall
Lx = 12.96                                                                                   # Streamwise length of the simulation
Lz = 2.2                                                                                     # Spanwise length of the simulation
symmetry = False                                                                             # Flag if the simulation can be considered symmetric in y-direction
delta_iteration = 5                                                                          # Iteration step for statistics computation

streamwise = 1                                                                               # Specify the streamwise direction used in MULTIFAST

preliminary_iterations = 150                                                                 # Number of preliminary iterations (by default set as start iteration)
number_iteration = 290                                                                       # Total number of iterations 
                                                                                             
################
custom_t = False                                                                            # Flag if the statistics are compute on a custom time range
rangeT = [130] 

root = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/'                             # Path that contains the raw results of the simulation
root_postproc = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/PostProcessing/'     # Path that will contain the post-processed results of the simulation

compute_lambda2 = False                                                                     # Flag if lambda2 will be computed
slice_lambda2 = 64                                                                          # Integer up to which is computed lambda2 (by default computing between 0 and slice_lambda2)

compute_spectrums = True                                                                    # Flag if spectrums will be computed  

################ For IBM_mask
IBM_flag = True 

# CHANNEL RE_TAU = 125, w*=12
i_start = [238,757,238,757,1276,1795,1276,1795] 
i_end = [280,799,280,799,1318,1837,1318,1837]    
j_start = [0,0,131,131,0,0,131,131] 
j_end = [48,48,180,180,48,48,180,180]   
k_start = [0,64,64,0,0,64,64,0] 
k_end = [63,128,128,63,63,128,128,63]   

x_start = [1.485,4.725,1.485,4.725,7.965,11.205,7.965,11.205] 
x_end = [1.755,4.995,1.755,4.995,8.235,11.475,8.235,11.475]    
y_start = [0,0,1.73,1.73,0,0,1.73,1.73] 
y_end = [0.27,0.27,2,2,0.27,0.27,2,2]   
z_start = [0,1.1,1.1,0,0,1.1,1.1,0] 
z_end = [1.1,2.2,2.2,1.1,1.1,2.2,2.2,1.1]   

###################################################################
################        CHOOSING THE SPOT           ###############
###################################################################

################ For turbulent spot
total_spot = False                                                                          # Flag if the total turbulent spot is used (mostly used in transitional flows)
heart_spot = False                                                                          # Flag if the core of the turbulent spot is used (mostly used in transitional flows)
wave_packet_spot = False                                                                    # Flag if the wave packet region of the turbulent spot is used (mostly used in transitional flows)

################ For any simulations
total_domain = True                                                                         # Flag if the entire domain is used as a spot

x_bounds = [15, 30]                                                             # Bounds of the spot in the streamwise direction
y_bounds = [0, 2]                                                                           # Bounds of the spot in the vertical direction
z_bounds = [0, 2.2]                                                                     # Bounds of the spot in the spanwise direction


###################################################################
################        STATS COMPUTATION           ###############
###################################################################
compute_velocity_stats = True                                                               # Flag if basic velocities stats are computed
compute_transport_equation_terms = True                                                     # Flag if budgets stats are computed
compute_reynolds_stresses = True                                                            # Flag if reynolds stresses stats are computed
compute_vorticity = True                                                                    # Flag if vorticity stats are computed
compute_temperature_stats = True                                                            # Flag if basic temperature stats are computed
compute_temperature_velocity_correlations = True                                            # Flag if temperature-velocity correlations stats are computed
compute_temperature_budgets = True                                                          # Flag if temperature budgets stats are computed

compute_E_tot_evolution = False                                                             # Flag if total kinetic energy is computed (mostly used in transitional flows)


###################################################################
################       CLASSIC PLOT INFOS           ###############
###################################################################

shift = 0                                                                                # Shift in the streamwise direction for better looking plot (mostly used in transitional flows)


################ 3D PLOTTING
plot_spot = False                                                                           # Flag if spot will be plotted (mostly used in transitional flows)
plot_thermal_spot = False                                                                   # Flag if thermal spot will be plotted (mostly used in transitional flows)
generate_video_spot = False                                                                 # Flag if videos of the spots are generated (mostly used in transitional flows)


################ 2D PLOTTING
plot_U_fluctuations = False                                                                 # Flag if fluctuations of U will be plotted
plot_V_fluctuations = False                                                                 # Flag if fluctuations of V will be plotted
plot_W_fluctuations = False                                                                 # Flag if fluctuations of W will be plotted
plot_T_fluctuations = False                                                                 # Flag if fluctuations of T will be plotted

plot_U_instantaneous = False                                                                # Flag if instantaneous of U will be plotted
plot_V_instantaneous = False                                                                # Flag if instantaneous of V will be plotted
plot_W_instantaneous = False                                                                # Flag if instantaneous of W will be plotted
plot_T_instantaneous = False                                                                # Flag if instantaneous of T will be plotted

plot_wall_shear_stress_contours = False                                                     # Flag if the wall shear stress will be plotted
plot_wall_temperature_flux_contours = False                                                 # Flag if the wall temperature flux will be plotted

plot_x_equal_to = False                                                                     # Location of the slice that wil be plotted in the ZY configuration (by default set as False)
plot_y_equal_to = False                                                                         # Location of the slice that wil be plotted in the XZ configuration (by default set as False)
plot_z_equal_to = 1.1                                                                     # Location of the slice that wil be plotted in the XY configuration (by default set as False)

colored_contours = False                                                                    # Flag if contours are colored with colormap
contours = True                                                                             # Flag if contours are plotted with lines


################ 1D PLOTTING
plot_velocity_stats = True                                                                  # Flag if basic velocities stats are plotted
plot_transport_equation_terms = True                                                        # Flag if budgets stats are plotted
plot_reynolds_stresses = True                                                               # Flag if reynolds stresses stats are plotted
plot_vorticity = True                                                                       # Flag if vorticity stats are plotted
plot_temperature_stats = True                                                               # Flag if basic temperature stats are plotted
plot_temperature_velocity_correlations = True                                               # Flag if temperature-velocity correlations stats are plotted
plot_temperature_budgets = True                                                             # Flag if temperature budgets stats are plotted

plot_E_tot_evolution = False                                                                # Flag if the temporal evolution of total kinetic energy is plotted (mostly used in transitional flows)
plot_Re_tau_evolution = False                                                               # Flag if the temporal evolution of Re_tau is plotted (mostly used in transitional flows)
plot_T_tau_evolution = False                                                                # Flag if the temporal evolution of T_tau is plotted (mostly used in transitional flows)
plot_Nu_evolution = False                                                                   # Flag if the temporal evolution of Nu is plotted (mostly used in transitional flows)


###################################################################
################      Henningson PLOT INFOS         ###############
##########    Check (DOI: 10.1017/S0022112093001429)      #########
###################################################################
plot_TEHL = False                                                                           # Flag if TEHL plots are used
