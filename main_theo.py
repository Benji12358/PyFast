"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Time-average statistics main

@description: Example of main for time-average statistics.
            Please, see the settings file for additional information.
"""

import sys
# appending a path
sys.path.append('./Utils')

from Settings import Settings
from CFD_mesh import CFD_mesh
from CFD_stats import CFD_stats

import plot_stats_utils

import warnings
warnings.filterwarnings("ignore")

compute_stats = 1
plot_stats = 1

def main():

    global settings, mesh, sim_environment    

    settings = Settings("postprocessSettings_Benj_theo")
    mesh = CFD_mesh(settings)
    
    sim_environment = CFD_stats(settings, mesh)
    
    sim_environment.get_spot(settings, mesh)
    if (settings.get_T):
        sim_environment.get_thermal_spot()
        
    if (compute_stats):
    
        sim_environment.compute_mean_fields(settings, mesh)
        
        print("Fields averaged!")
            
        sim_environment.init_ui_stats_arrays(settings, mesh)
        if (settings.get_T):
            sim_environment.init_T_stats_arrays(settings, mesh)
        
        print("Computing statistics ...")
        
        for t in settings.all_iterations:
            
            sim_environment.read_current_fields(settings, t)
            
            ###################
            # for dynamic stats
            ###################
            sim_environment.compute_basic_stats()
            sim_environment.compute_vorticity(settings, mesh)
            sim_environment.compute_transport_equation_terms(settings, mesh)
            
            if (settings.get_T):
                sim_environment.compute_basic_thermal_stats()
                sim_environment.compute_transport_equation_terms_temperature(settings, mesh)
        
        print("Statistics computed!")
    
        sim_environment.average_ui_stats_arrays(settings)
        if (settings.get_T):
            sim_environment.average_T_stats_arrays(settings)
            
        sim_environment.save_basic_velocity_stats(settings, mesh)
        sim_environment.save_reynolds_stresses(settings, mesh)
        sim_environment.save_vorticity(settings, mesh)
        sim_environment.save_transport_equation_terms(settings, mesh)
        
        if (settings.get_T):
            sim_environment.save_basic_temperature_stats(settings, mesh)
            sim_environment.save_transport_equation_terms_temperature(settings, mesh)
        
        print("Statistics saved!")


	#plot stats
    
    if (plot_stats):
    
        plot_stats_utils.read_and_plot_basic_velocity_stats(settings, mesh, plotWithTheo=False)
        plot_stats_utils.read_and_plot_transport_equation_terms(settings, mesh, plotWithTheo=False)
        plot_stats_utils.read_and_plot_vorticity(settings, mesh, plotWithTheo=False)
        plot_stats_utils.read_and_plot_reynolds_stresses(settings, mesh, plotWithTheo=False)
        
        plot_stats_utils.read_and_plot_basic_temperature_stats(settings, mesh, plotWithTheo=False)
        plot_stats_utils.read_and_plot_temperature_budgets(settings, mesh, plotWithTheo=False)


	#plots fields





main()