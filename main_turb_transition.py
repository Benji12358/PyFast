"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Turbulent spot statistics main

@description: Example of main for statistics in a turbulent spot.
            Please, see the settings file for additional information.
"""

import sys
# appending a path
sys.path.append('./Utils')

from Settings import Settings
from CFD_mesh import CFD_mesh
from CFD_stats import CFD_stats

from plot_utils import plot_2D_contours_with_spot

import plot_stats_utils

import warnings
warnings.filterwarnings("ignore")

def main():

    global settings, mesh, sim_environment    

    """Read the settings define by the user, and build the mesh"""
    settings = Settings("postprocessSettings_Benj")
    mesh = CFD_mesh(settings)
    
    """Initialize all variables of the simulation"""
    sim_environment = CFD_stats(settings, mesh)

    """Initialize the mean fields"""
    sim_environment.set_mean_fields(settings, mesh)

    if (settings.plot_TE):
        
        for t in settings.all_iterations:
            
            """Read the instantaneous fields"""
            sim_environment.read_current_fields(settings, t)
            
            """Initialize statistics arrays"""
            sim_environment.init_ui_stats_arrays(settings, mesh)
            if (settings.get_T):
                sim_environment.init_T_stats_arrays(settings, mesh)
            
            if (settings.get_dynamic_spot):
                sim_environment.get_turbulent_spot(settings, mesh)
                
            print("Computing statistics ...")
            
            """Compute and save velocities statistics"""
            if (settings.get_dynamic_stats):
            
                sim_environment.set_mean_spot(settings, mesh)
                sim_environment.compute_basic_stats()
                    
                if (settings.compute_velocity_stats):
                    sim_environment.save_basic_velocity_stats(settings, mesh, t)
                    
                if (settings.compute_reynolds_stresses):
                    sim_environment.save_reynolds_stresses(settings, mesh, t)
                    
                if (settings.compute_vorticity):
                    sim_environment.compute_vorticity(settings, mesh)
                    sim_environment.save_vorticity(settings, mesh, t)
                
                if (settings.compute_transport_equation_terms):
                    sim_environment.compute_transport_equation_terms(settings, mesh)
                    sim_environment.save_transport_equation_terms(settings, mesh, t)
                    
            """Compute and save Lambda2"""
            if (settings.compute_lambda2): 
                sim_environment.compute_lambda2(mesh, settings, sim_environment.Re_tau_spot, t)
            
            """Compute and save temperature statistics"""
            if (settings.get_thermal_spot):
                sim_environment.get_turbulent_thermal_spot(settings, mesh)
            
            if (settings.get_thermal_stats):
            
                sim_environment.set_mean_thermal_spot(settings, mesh)
                sim_environment.compute_basic_thermal_stats()
                    
                if (settings.compute_temperature_stats):
                    sim_environment.save_basic_temperature_stats(settings, mesh, t)
                
                if (settings.compute_temperature_budgets):
                    sim_environment.compute_transport_equation_terms_temperature(settings, mesh)
                    sim_environment.save_transport_equation_terms_temperature(settings, mesh, t)
    
            """Plot the velocities statistics"""
            plot_stats_utils.read_and_plot_basic_velocity_stats(settings, mesh, t, plotWithTheo=True)
            plot_stats_utils.read_and_plot_transport_equation_terms(settings, mesh, t, plotWithTheo=True)
            plot_stats_utils.read_and_plot_vorticity(settings, mesh, t, plotWithTheo=True)
            plot_stats_utils.read_and_plot_reynolds_stresses(settings, mesh, t, plotWithTheo=True)
            
            """Plot the temperature statistics"""
            plot_stats_utils.read_and_plot_basic_temperature_stats(settings, mesh, t, plotWithTheo=True)
            plot_stats_utils.read_and_plot_temperature_budgets(settings, mesh, t, plotWithTheo=True)
    
            """Plot wall shear stress contours"""
            if (settings.plot_wall_shear_stress_contours):
                
                wall_ST_fluc = sim_environment.get_wall_shear_stress(mesh)
                plot_2D_contours_with_spot(wall_ST_fluc, sim_environment.spot[:,0,:], mesh, settings)
                
            """Plot wall temperature flux contours"""
            if (settings.plot_wall_temperature_flux_contours):
                
                wall_TF_fluc = sim_environment.get_wall_temperature_flux(mesh)
                plot_2D_contours_with_spot(wall_TF_fluc, sim_environment.thermal_spot[:,0,:], mesh, settings)
        
            """Plot 2D fields"""
            sim_environment.plot_2D_fields(settings, mesh)



main()