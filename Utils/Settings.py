"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class Settings

@description: Read a settings file written in Python and gather all flags/variables.
            The object contains all settings of the user.
"""

import numpy as np

class Settings:
    
    def __init__(self, path_settings):
        """Initialisation of the Settings object.
        
        Parameters
        ----------
        path_settings : String
            Path of the settings file to read, to get all user settings
    
        Returns
        -------
        Settings
            An object containing all user settings.
        """

        # read settings
        exec(open(path_settings + ".py").read(), globals())
        
        ###################################################################
        ################    GENERAL SIMULATION INFOS        ###############
        ###################################################################
        self.get_main_infos()

        ###################################################################
        ################        CHOOSING THE SPOT           ###############
        ###################################################################
        self.get_spot_infos()

        ###################################################################
        ###################### Classical PLOT INFOS #######################
        ###################################################################
        self.shift = shift
        self.get_stats_computation_infos()

        self.get_2D_infos()
        self.get_3D_infos()
        self.get_stats_plotting_infos()
              
        self.plot_TEHL = plot_TEHL
        
        if (self.plot_TEHL):
            self.get_TEHL_info()

        # compile all information
        self.gather_information()
        self.get_time_range()


    def gather_information(self):
        """Initialisation of key booleans used in the main environment.
    
            Update the key boolean of a Settings object.
        """

        self.get_dynamic_stats = self.compute_velocity_stats or self.compute_transport_equation_terms \
                            or self.compute_reynolds_stresses or self.compute_vorticity

        self.get_thermal_stats = self.compute_temperature_stats or self.compute_temperature_velocity_correlations \
                            or self.compute_temperature_budgets

        self.plot_stats = self.plot_velocity_stats or self.plot_transport_equation_terms \
                            or self.plot_reynolds_stresses or self.plot_vorticity \
                            or self.plot_temperature_stats or self.plot_temperature_velocity_correlations \
                            or self.plot_temperature_budgets

        self.get_T = self.plot_T_instantaneous or self.plot_T_fluctuations \
                            or self.plot_wall_temperature_flux_contours \
                            or self.compute_temperature_stats \
                            or self.compute_temperature_velocity_correlations or self.plot_temperature_stats \
                            or self.plot_temperature_velocity_correlations or self.plot_thermal_spot

        self.get_dynamic_spot = self.plot_E_tot_evolution or self.plot_spot \
                            or self.plot_wall_shear_stress_contours or self.get_dynamic_stats

        self.get_thermal_spot = self.plot_thermal_spot or self.get_thermal_stats \
                            or self.plot_wall_temperature_flux_contours

        self.plot_TE = self.plot_U_fluctuations or self.plot_V_fluctuations \
                            or self.plot_W_fluctuations or self.plot_U_instantaneous \
                            or self.plot_V_instantaneous or self.plot_W_instantaneous \
                            or self.plot_T_fluctuations or self.plot_T_instantaneous \
                            or self.plot_spot or self.plot_thermal_spot \
                            or self.plot_wall_shear_stress_contours or self.plot_wall_temperature_flux_contours \
                            or self.compute_lambda2 or self.get_dynamic_stats \
                            or self.get_thermal_stats or self.plot_stats


    def get_main_infos(self):
        """Initialization of the main information of the simulation to post-process
        """

        # initialize main parameters for the simulation
        self.current_path = current_path
        self.mean_path = mean_path
        self.theo_path = theo_path
        self.Re = Re
        self.Pr = Pr
        self.deltaT = deltaT
        self.Lx = Lx
        self.Lz = Lz
        self.symmetry = symmetry
        self.delta_iteration = delta_iteration

        # initialize number of iterations
        self.preliminary_iterations = preliminary_iterations
        self.number_iteration = number_iteration

        self.custom_t = custom_t
        self.rangeT = rangeT

        # initialize path for post-processing
        self.root = root
        self.root_postproc = root_postproc
        self.path_postProcess = self.root_postproc + self.current_path

        # initialize boolean for lambda2 computing
        self.compute_lambda2 = compute_lambda2
        self.slice_lambda2 = slice_lambda2


    def get_spot_infos(self):
        """Initialization of the main information of the spot
        """
        
        # initialize parameters for spot/domain in which stats are computed
        self.total_spot = total_spot
        self.heart_spot = heart_spot
        self.wave_packet_spot = wave_packet_spot

        self.total_domain = total_domain
        
        self.x_bounds = x_bounds
        self.y_bounds = y_bounds
        self.z_bounds = z_bounds


    def get_stats_computation_infos(self):
        """Initialization of the main information for stats computation
        """

        self.compute_velocity_stats = compute_velocity_stats
        self.compute_transport_equation_terms = compute_transport_equation_terms
        self.compute_reynolds_stresses = compute_reynolds_stresses
        self.compute_vorticity = compute_vorticity

        self.compute_temperature_stats = compute_temperature_stats
        self.compute_temperature_velocity_correlations = compute_temperature_velocity_correlations
        self.compute_temperature_budgets = compute_temperature_budgets

        self.compute_E_tot_evolution = compute_E_tot_evolution


    def get_stats_plotting_infos(self):
        """Initialization of the main information for stats plotting
        """

        self.plot_velocity_stats = plot_velocity_stats
        self.plot_transport_equation_terms = plot_transport_equation_terms
        self.plot_reynolds_stresses = plot_reynolds_stresses
        self.plot_vorticity = plot_vorticity

        self.plot_temperature_stats = plot_temperature_stats
        self.plot_temperature_velocity_correlations = plot_temperature_velocity_correlations
        self.plot_temperature_budgets = plot_temperature_budgets

        self.plot_E_tot_evolution = plot_E_tot_evolution
        self.plot_Re_tau_evolution = plot_Re_tau_evolution
        self.plot_T_tau_evolution = plot_T_tau_evolution
        self.plot_Nu_evolution = plot_Nu_evolution


    def get_2D_infos(self):
        """Initialization of the main information for 2D plotting
        """

        self.plot_U_fluctuations = plot_U_fluctuations
        self.plot_V_fluctuations = plot_V_fluctuations
        self.plot_W_fluctuations = plot_W_fluctuations
        self.plot_T_fluctuations = plot_T_fluctuations

        self.plot_U_instantaneous = plot_U_instantaneous
        self.plot_V_instantaneous = plot_V_instantaneous
        self.plot_W_instantaneous = plot_W_instantaneous
        self.plot_T_instantaneous = plot_T_instantaneous

        self.plot_wall_shear_stress_contours = plot_wall_shear_stress_contours
        self.plot_wall_temperature_flux_contours = plot_wall_temperature_flux_contours

        self.plot_x_equal_to = plot_x_equal_to
        self.plot_y_equal_to = plot_y_equal_to
        self.plot_z_equal_to = plot_z_equal_to
        
        self.colored_contours = colored_contours
        self.contours = contours


    def get_3D_infos(self):
        """Initialization of the main information for 3D plotting
        """

        self.plot_spot = plot_spot
        self.plot_thermal_spot = plot_thermal_spot
        self.generate_video_spot = generate_video_spot


    def get_time_range(self):
        """Initialisation of the array containing all timestep of the simulation to post process.
    
            Update the iteration array of a Settings object.
        """

        if (self.custom_t):
            self.all_iterations = [x + self.preliminary_iterations for x in self.rangeT]
        else:
            self.all_iterations = range(self.preliminary_iterations,self.number_iteration+1,self.delta_iteration)
            
        self.all_iterations = np.array(self.all_iterations)


    def get_TEHL_info(self):
        """Get all information related to TEHL plots
           This is to plot temporal evolution like Henningson Fig. 2 (DOI: 10.1017/S0022112093001429)
           Mostly used for transitionnal flows
        """

        # initialize boolean for Henningson like 2D plotting
        self.plot_TEHL_U_fluctuations = plot_TEHL_U_fluctuations
        self.plot_TEHL_V_fluctuations = plot_TEHL_V_fluctuations
        self.plot_TEHL_W_fluctuations = plot_TEHL_W_fluctuations
        self.plot_TEHL_T_fluctuations = plot_TEHL_T_fluctuations

        self.plot_TEHL_wall_shear_stress_contours = plot_TEHL_wall_shear_stress_contours
        self.plot_TEHL_wall_temperature_flux_contours = plot_TEHL_wall_temperature_flux_contours

        # initialize boolean for Henningson like statistics plotting
        self.plot_TEHL_velocity_stats = plot_TEHL_velocity_stats
        self.plot_TEHL_transport_equation_terms = plot_TEHL_transport_equation_terms
        self.plot_TEHL_reynolds_stresses = plot_TEHL_reynolds_stresses
        
        self.plot_TEHL_temperature_stats = plot_TEHL_temperature_stats
        self.plot_TEHL_temperature_velocity_correlations = plot_TEHL_temperature_velocity_correlations

        # initialize data for Henningson like plots
        self.time_values = time_values
        self.y_limits = y_limits
        self.x_limits = x_limits