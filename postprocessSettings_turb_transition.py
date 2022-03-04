###################################################################
################    GENERAL SIMULATION INFOS        ###############
###################################################################

# L=32*Pi*H / epsilon=0.14 / theta=0°
current_path = 'PL_Vortices_LargeEps_fine'                    								# Name of the simulation to post-process
mean_path = 'Preliminary_Laminar_PFlow'														# Name of the simulation containing mean fields (by default, put same name as simulation to post-process)
theo_path = 'ReTau_150_PS_Dir'																# Name of the simulation containing theoretical results (by default, leave blank)
Re = 3000																					# Reynolds number of the simulation
Pr = 1																						# Prandtl number of the simulation
deltaT = 0.5                                                                                # Temperature step between the lower and the upper wall
Lx = 32*np.pi                                                                               # Streamwise length of the simulation
Lz = 8*np.pi                                                                                # Spanwise length of the simulation
symmetry = True                                                                             # Flag if the simulation can be considered symmetric in y-direction
delta_iteration = 5                                                                         # Iteration step for statistics computation

streamwise = 1                                                                         		# Specify the streamwise direction used in MULTIFAST

preliminary_iterations = 10                                                                 # Number of preliminary iterations (by default set as start iteration)
number_iteration = 175                                                                      # Total number of iterations 

custom_t = True                                                                            	# Flag if the statistics are compute on a custom time range
rangeT = [165]                                                                              # Custom time range array (by default set anything)

root = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/'                             # Path that contains the raw results of the simulation
root_postproc = '/fsnet/people/arrondea7b/project/21MULTIFAST/Sim_data/PostProcessing/'     # Path that will contain the post-processed results of the simulation

compute_lambda2 = False                                                                     # Flag if lambda2 will be computed
slice_lambda2 = 64                                                                          # Integer up to which is computed lambda2 (by default computing between 0 and slice_lambda2)


###################################################################
################        CHOOSING THE SPOT           ###############
###################################################################

################ For turbulent spot
total_spot = True                                                                          	# Flag if the total turbulent spot is used (mostly used in transitional flows)
heart_spot = False                                                                          # Flag if the core of the turbulent spot is used (mostly used in transitional flows)
wave_packet_spot = False                                                                    # Flag if the wave packet region of the turbulent spot is used (mostly used in transitional flows)

################ For any simulations
total_domain = False                                                                        # Flag if the entire domain is used as a spot

x_bounds = [24*np.pi, 32*np.pi]                                                             # Bounds of the spot in the streamwise direction
y_bounds = [0, 2]                                                                           # Bounds of the spot in the vertical direction
z_bounds = [0, 8*np.pi]                                                                     # Bounds of the spot in the spanwise direction


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

shift = 77.5                                                                                # Shift in the streamwise direction for better looking plot (mostly used in transitional flows)


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

plot_wall_shear_stress_contours = True                                                     	# Flag if the wall shear stress will be plotted
plot_wall_temperature_flux_contours = True                                                 	# Flag if the wall temperature flux will be plotted

plot_x_equal_to = False                                                                     # Location of the slice that wil be plotted in the ZY configuration (by default set as False)
plot_y_equal_to = 1                                                                         # Location of the slice that wil be plotted in the XZ configuration (by default set as False)
plot_z_equal_to = False                                                                     # Location of the slice that wil be plotted in the XY configuration (by default set as False)

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

# plot temporal evolution like Henningson Fig 2
plot_TEHL_U_fluctuations = False															# Flag if fluctuations of U will be plotted in TEHL way
plot_TEHL_V_fluctuations = False															# Flag if fluctuations of V will be plotted in TEHL way
plot_TEHL_W_fluctuations = False															# Flag if fluctuations of W will be plotted in TEHL way
plot_TEHL_T_fluctuations = False															# Flag if fluctuations of T will be plotted in TEHL way

plot_TEHL_wall_shear_stress_contours = False												# Flag if the wall shear stress will be plotted in TEHL way
plot_TEHL_wall_temperature_flux_contours = False											# Flag if the wall temperature flux will be plotted in TEHL way

# 1D stats plotting
plot_TEHL_velocity_stats = False															# Flag if basic velocities stats are plotted in TEHL way
plot_TEHL_transport_equation_terms = False													# Flag if budgets are plotted in TEHL way
plot_TEHL_reynolds_stresses = False															# Flag if reynolds stresses stats are plotted in TEHL way
plot_TEHL_temperature_stats = False															# Flag if basic temperature stats are plotted in TEHL way
plot_TEHL_temperature_velocity_correlations = False											# Flag if temperature-velocity correlations stats are plotted in TEHL way

time_values = [89,109,129,149]																# Time values at which the stats are plotted
y_limits = [[-15,15], [-15,15], [-15,15], [-15,15]]											# Y-axis bounds for 2D plots in TEHL way
x_limits = [[10,45], [15,50], [20,60], [25,70]]												# X-axis bounds for 2D plots in TEHL way