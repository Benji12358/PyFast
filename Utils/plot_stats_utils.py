"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Utilities functions for 1D plotting

@description: Contains general functions for 1D plotting for classical statistical data 
            (Mean, RMS, Skewness, Flatness, Vorticity, Reynolds stresses, Transport equation terms) for the velocities
            (Mean, RMS, Skewness, Flatness, Correlations, Transport equation terms) for the temperature
"""

import numpy as np

from CFD_plot import CFD_plot
from CFD_mesh import CFD_mesh
from Settings import Settings


from utilities import find_nearest, print_infos
from stats_utils import read_stats, read_U_main_quantities, read_T_main_quantities 
        
def read_and_plot_basic_velocity_stats(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the mean, RMS, skewness and flatness velocities. Plot the diagnostic function.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_postProcess = settings.path_postProcess+'/VelocityStats/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    
    wall_UQ = read_U_main_quantities(current_path_postProcess, current_iteration)
        
    # global Mean_RMS
    Mean_RMS = {
        'name': 'Mean_RMS',
        
        'Wall normal coordinate': None,
        'U_mean': None,
        'V_mean': None,
        'W_mean': None,
        'U_rms': None,
        'V_rms': None,
        'W_rms': None
        }
    
    Mean_RMS = read_stats(Mean_RMS, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_UQ_theo = read_U_main_quantities(current_path_postProcess_theo)
    
        # global Mean_RMS
        Mean_RMS_theo = {
            'name': 'Mean_RMS',
            
            'Wall normal coordinate': None,
            'U_mean': None,
            'V_mean': None,
            'W_mean': None,
            'U_rms': None,
            'V_rms': None,
            'W_rms': None
            }
        
        Mean_RMS_theo = read_stats(Mean_RMS_theo, current_path_postProcess_theo)
    
    ######## mean
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_mean']*(settings.Re/wall_UQ['Re_tau_spot']), 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_mean']*(settings.Re/wall_UQ['Re_tau_spot']), linestyle=':', linewidth=2.5, color='k')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_RMS_theo['U_mean']*(settings.Re/wall_UQ_theo['Re_tau_spot']))
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_mean']*(settings.Re/wall_UQ['Re_tau_spot']))
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'],axis_scale='log')
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=5*10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'],axis_scale='log')
    fig.chg_y_axis(r'$\overline{U}^{+}$',axis_low_bound=0,axis_high_bound=20,axis_ticks = [0,5,10,15,20])
    fig.custom_layout()
    fig.add_ticks_x(axis_scale='log')
    fig.add_ticks_y()
    fig.display()

    ####### RMS
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']*(settings.Re/wall_UQ['Re_tau_spot']), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$\overline{u^{'+}_{RMS}}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['V_rms']*(settings.Re/wall_UQ['Re_tau_spot']), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$\overline{v^{'+}_{RMS}}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['W_rms']*(settings.Re/wall_UQ['Re_tau_spot']), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\overline{w^{'+}_{RMS}}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']*(settings.Re/wall_UQ['Re_tau_spot']), linestyle=':', linewidth=2.5, color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['V_rms']*(settings.Re/wall_UQ['Re_tau_spot']), linestyle=':', linewidth=2.5, color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['W_rms']*(settings.Re/wall_UQ['Re_tau_spot']), linestyle=':', linewidth=2.5, color='tab:green')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_RMS_theo['U_rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']), color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_RMS_theo['V_rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']), color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_RMS_theo['W_rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']), color='tab:green')
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:blue', label_name=r"$\overline{u^{'+}_{RMS}}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['V_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:orange', label_name=r"$\overline{v^{'+}_{RMS}}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['W_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:green', label_name=r"$\overline{w^{'+}_{RMS}}$")
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=2*wall_UQ['Re_tau_spot'])
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=5*10**(-3),axis_high_bound=10**0,axis_scale='log')
    fig.chg_y_axis('Velocities RMS',axis_low_bound=0,axis_high_bound=1.1*np.max(Mean_RMS['U_rms']*settings.Re/wall_UQ['Re_tau_spot']))
    fig.custom_layout(enableLegend=True)
    fig.add_ticks_x()
    fig.add_ticks_y()
    fig.display()
    
    # fig = CFD_plot('full')
    # fig.add_plot(mesh.Yc, Mean_RMS['U_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:blue', label_name=r"$\overline{u^{'+}_{RMS}}$")
    # fig.add_plot(mesh.Yc, Mean_RMS['V_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:orange', label_name=r"$\overline{v^{'+}_{RMS}}$")
    # fig.add_plot(mesh.Yc, Mean_RMS['W_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:green', label_name=r"$\overline{w^{'+}_{RMS}}$")
    # # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_x_axis(r'$y^{*}$',axis_low_bound=5*10**(-3),axis_high_bound=10**0,axis_scale='log')
    # fig.chg_y_axis('Velocities RMS',axis_low_bound=0,axis_high_bound=1.1*np.max(Mean_RMS['U_rms']*settings.Re/wall_UQ['Re_tau_spot']))
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    # fig = CFD_plot('full')
    # fig.add_plot(mesh.Yc, Mean_RMS['U_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:blue', label_name=r"$\overline{u^{'+}_{RMS}}$")
    # fig.add_plot(mesh.Yc, Mean_RMS['V_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:orange', label_name=r"$\overline{v^{'+}_{RMS}}$")
    # fig.add_plot(mesh.Yc, Mean_RMS['W_rms']*(settings.Re/wall_UQ['Re_tau_spot']), color='tab:green', label_name=r"$\overline{w^{'+}_{RMS}}$")
    # # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_x_axis(r'$y^{*}$',axis_low_bound=0,axis_high_bound=1)
    # fig.chg_y_axis('Velocities RMS',axis_low_bound=0,axis_high_bound=1.1*np.max(Mean_RMS['U_rms']*settings.Re/wall_UQ['Re_tau_spot']))
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()

    # ######## diagnostic functions
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']/Mean_RMS['U_mean'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$\frac{\overline{u^{'+}_{RMS}}}{\overline{u}}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']/Mean_RMS['U_mean'], linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_RMS_theo['U_rms']/Mean_RMS_theo['U_mean'], color='tab:blue')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_RMS['U_rms']/Mean_RMS['U_mean'], color='tab:blue', label_name=r"$\frac{\overline{u^{'+}_{RMS}}}{\overline{u}}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_y_axis('')
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
    
    # Skewness_Flatness = {
    #     'name': 'Skewness_Flatness',
        
    #     'Wall normal coordinate': None,
    #     'U_skewness': None,
    #     'V_skewness': None,
    #     'W_skewness': None,
    #     'U_flatness': None,
    #     'V_flatness': None,
    #     'W_flatness': None
    #     }
    
    # Skewness_Flatness = read_stats(Skewness_Flatness, current_path_postProcess, current_iteration) 
    
    # if (plotWithTheo):
    
    #     # global Mean_RMS
    #     Skewness_Flatness_theo = {
    #         'name': 'Skewness_Flatness',
            
    #         'Wall normal coordinate': None,
    #         'U_skewness': None,
    #         'V_skewness': None,
    #         'W_skewness': None,
    #         'U_flatness': None,
    #         'V_flatness': None,
    #         'W_flatness': None
    #         }
        
    #     Skewness_Flatness_theo = read_stats(Skewness_Flatness_theo, current_path_postProcess_theo)

    # ######## Skewness
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_skewness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$S(u^{'+})$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_skewness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$S(v^{'+})$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_skewness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$S(w^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_skewness'], linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_skewness'], linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_skewness'], linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['U_skewness'], color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['V_skewness'], color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['W_skewness'], color='tab:green')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_skewness'], color='tab:blue', label_name=r"$S(u^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_skewness'], color='tab:orange', label_name=r"$S(v^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_skewness'], color='tab:green', label_name=r"$S(w^{'+})$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_y_axis('Velocities skewness')
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()

    # ######## Flatness
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_flatness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$F(u^{'+})$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_flatness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$F(v^{'+})$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_flatness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$F(w^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_flatness'], linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_flatness'], linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_flatness'], linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['U_flatness'], color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['V_flatness'], color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['W_flatness'], color='tab:green')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['U_flatness'], color='tab:blue', label_name=r"$F(u^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['V_flatness'], color='tab:orange', label_name=r"$F(v^{'+})$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['W_flatness'], color='tab:green', label_name=r"$F(w^{'+})$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_y_axis('Velocities flatness')
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
        
def read_and_plot_transport_equation_terms(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the budgets of uu, uv, vv, ww and k.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_vel_stats_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/TransportEquationsTerms/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/TransportEquationsTerms/'
    
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    ####################################################
    #dUfUf/dt
    uu_budget = {
        'name': 'uu_budget',
        'Wall normal coordinate': None,
        'Production': None,
        'Turbulent Diffusive Flux': None,
        'Pressure velocity gradient correlation': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    uu_budget = read_stats(uu_budget, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_UQ_theo = read_U_main_quantities(current_path_vel_stats_theo)
    
        # global Mean_RMS
        uu_budget_theo = {
            'name': 'uu_budget',
            
            'Wall normal coordinate': None,
            'Production': None,
            'Turbulent Diffusive Flux': None,
            'Pressure velocity gradient correlation': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
        
        uu_budget_theo = read_stats(uu_budget_theo, current_path_postProcess_theo)
    
    # # plot
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:purple')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uu_budget_theo['Production']                             * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uu_budget_theo['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uu_budget_theo['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uu_budget_theo['Molecular diffusion']                    * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uu_budget_theo['Dissipation']                            * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:purple')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:red', label_name=r"$D^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uu_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:purple', label_name=r"$\epsilon^{+}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    # fig.chg_y_axis('')
    # fig.custom_layout(enableLegend=True)
    # fig.add_title(r'$uu$' + ' budget')
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    ####################################################
    #dVfVf/dt
    vv_budget = {
        'name': 'vv_budget',
        'Wall normal coordinate': None,
        'Production': None,
        'Turbulent Diffusive Flux': None,
        'Pressure velocity gradient correlation': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    vv_budget = read_stats(vv_budget, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        # global Mean_RMS
        vv_budget_theo = {
            'name': 'vv_budget',
            
            'Wall normal coordinate': None,
            'Production': None,
            'Turbulent Diffusive Flux': None,
            'Pressure velocity gradient correlation': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
        
        vv_budget_theo = read_stats(vv_budget_theo, current_path_postProcess_theo)
    
    # # plot
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:purple')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vv_budget_theo['Production']                             * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vv_budget_theo['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vv_budget_theo['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vv_budget_theo['Molecular diffusion']                    * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vv_budget_theo['Dissipation']                            * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:purple')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:red', label_name=r"$D^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:purple', label_name=r"$\epsilon^{+}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    # fig.chg_y_axis('')
    # fig.custom_layout(enableLegend=True)
    # fig.add_title(r'$vv$' + ' budget')
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    ####################################################
    #dWfWf/dt
    ww_budget = {
        'name': 'ww_budget',
        'Wall normal coordinate': None,
        'Production': None,
        'Turbulent Diffusive Flux': None,
        'Pressure velocity gradient correlation': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    ww_budget = read_stats(ww_budget, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        # global Mean_RMS
        ww_budget_theo = {
            'name': 'ww_budget',
            
            'Wall normal coordinate': None,
            'Production': None,
            'Turbulent Diffusive Flux': None,
            'Pressure velocity gradient correlation': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
        
        ww_budget_theo = read_stats(ww_budget_theo, current_path_postProcess_theo)
    
    # # plot
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:purple')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], ww_budget_theo['Production']                             * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], ww_budget_theo['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], ww_budget_theo['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], ww_budget_theo['Molecular diffusion']                    * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], ww_budget_theo['Dissipation']                            * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:purple')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:red', label_name=r"$D^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], ww_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:purple', label_name=r"$\epsilon^{+}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    # fig.chg_y_axis('')
    # fig.custom_layout(enableLegend=True)
    # fig.add_title(r'$ww$' + ' budget')
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    ####################################################
    #dUfVf/dt
    uv_budget = {
        'name': 'uv_budget',
        'Wall normal coordinate': None,
        'Production': None,
        'Turbulent Diffusive Flux': None,
        'Pressure velocity gradient correlation': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    uv_budget = read_stats(uv_budget, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        # global Mean_RMS
        uv_budget_theo = {
            'name': 'uv_budget',
            
            'Wall normal coordinate': None,
            'Production': None,
            'Turbulent Diffusive Flux': None,
            'Pressure velocity gradient correlation': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
        
        uv_budget_theo = read_stats(uv_budget_theo, current_path_postProcess_theo)
    
    # # plot
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:purple')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uv_budget_theo['Production']                             * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:blue')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uv_budget_theo['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:orange')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uv_budget_theo['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:green')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uv_budget_theo['Molecular diffusion']                    * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:red')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], uv_budget_theo['Dissipation']                            * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:purple')
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:blue', label_name=r"$P^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:orange', label_name=r"$T^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:green', label_name=r"$\pi^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:red', label_name=r"$D^{+}$")
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], uv_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:purple', label_name=r"$\epsilon^{+}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    # fig.chg_y_axis('')
    # fig.custom_layout(enableLegend=True)
    # fig.add_title(r'$uv$' + ' budget')
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    ####################################################
    #dk/dt
    k_budget = {
        'name': 'k_budget',
        'Wall normal coordinate': None,
        'Production': None,
        'Turbulent Diffusive Flux': None,
        'Pressure velocity gradient correlation': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    k_budget['Production'] = (uu_budget['Production'] + vv_budget['Production'] + ww_budget['Production'])/2
    k_budget['Turbulent Diffusive Flux'] = (uu_budget['Turbulent Diffusive Flux'] + vv_budget['Turbulent Diffusive Flux'] + ww_budget['Turbulent Diffusive Flux'])/2
    k_budget['Pressure velocity gradient correlation'] = (uu_budget['Pressure velocity gradient correlation'] + vv_budget['Pressure velocity gradient correlation'] + ww_budget['Pressure velocity gradient correlation'])/2
    k_budget['Molecular diffusion'] = (uu_budget['Molecular diffusion'] + vv_budget['Molecular diffusion'] + ww_budget['Molecular diffusion'])/2
    k_budget['Dissipation'] = (uu_budget['Dissipation'] + vv_budget['Dissipation'] + ww_budget['Dissipation'])/2
    
    if (plotWithTheo):
    
        # global Mean_RMS
        k_budget_theo = {
            'name': 'k_budget',
            
            'Wall normal coordinate': None,
            'Production': None,
            'Turbulent Diffusive Flux': None,
            'Pressure velocity gradient correlation': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
    
        k_budget_theo['Production'] = (uu_budget_theo['Production'] + vv_budget_theo['Production'] + ww_budget_theo['Production'])/2
        k_budget_theo['Turbulent Diffusive Flux'] = (uu_budget_theo['Turbulent Diffusive Flux'] + vv_budget_theo['Turbulent Diffusive Flux'] + ww_budget_theo['Turbulent Diffusive Flux'])/2
        k_budget_theo['Pressure velocity gradient correlation'] = (uu_budget_theo['Pressure velocity gradient correlation'] + vv_budget_theo['Pressure velocity gradient correlation'] + ww_budget_theo['Pressure velocity gradient correlation'])/2
        k_budget_theo['Molecular diffusion'] = (uu_budget_theo['Molecular diffusion'] + vv_budget_theo['Molecular diffusion'] + ww_budget_theo['Molecular diffusion'])/2
        k_budget_theo['Dissipation'] = (uu_budget_theo['Dissipation'] + vv_budget_theo['Dissipation'] + ww_budget_theo['Dissipation'])/2
        
    # plot
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\pi^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:green')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:red')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), linestyle=':', linewidth=2.5, color='tab:purple')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], k_budget_theo['Production']                             * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], k_budget_theo['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], k_budget_theo['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:green')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], k_budget_theo['Molecular diffusion']                    * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:red')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], k_budget_theo['Dissipation']                            * (settings.Re**3/wall_UQ_theo['Re_tau_spot']**4), color='tab:purple')
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Production']                             * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:blue', label_name=r"$P^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Turbulent Diffusive Flux']               * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:orange', label_name=r"$T^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Pressure velocity gradient correlation'] * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:green', label_name=r"$\pi^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Molecular diffusion']                    * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:red', label_name=r"$D^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], k_budget['Dissipation']                            * (settings.Re**3/wall_UQ['Re_tau_spot']**4), color='tab:purple', label_name=r"$\epsilon^{+}$")
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    fig.chg_y_axis('')
    fig.custom_layout(enableLegend=True)
    fig.add_title(r'$k$' + ' budget')
    fig.add_ticks_x(axis_scale='log')
    fig.add_ticks_y()
    fig.display()
    
    
def read_and_plot_vorticity(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the rms of the vorticity.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_vel_stats_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/Vorticity/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/Vorticity/'
    
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    ####################################################
    #dUfUf/dt
    vorticity = {
        'name': 'Vort_fluc',
            
        'Wall normal coordinate': None,
        'omega_x': None,
        'omega_y': None,
        'omega_z': None,
        'omega_x rms': None,
        'omega_y rms': None,
        'omega_z rms': None
        }
    
    vorticity = read_stats(vorticity, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_UQ_theo = read_U_main_quantities(current_path_vel_stats_theo)
    
        # global Mean_RMS
        vorticity_theo = {
            'name': 'Vort_fluc',
                
            'Wall normal coordinate': None,
            'omega_x': None,
            'omega_y': None,
            'omega_z': None,
            'omega_x rms': None,
            'omega_y rms': None,
            'omega_z rms': None
            }
        
        vorticity_theo = read_stats(vorticity_theo, current_path_postProcess_theo)
    
    # vorticity RMS
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_x rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$\omega_{x~RMS}^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_y rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$\omega_{y~RMS}^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_z rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\omega_{z~RMS}^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_x rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), linestyle=':', linewidth=2.5, color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_y rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), linestyle=':', linewidth=2.5, color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_z rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), linestyle=':', linewidth=2.5, color='tab:green')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vorticity_theo['omega_x rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']**2), color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vorticity_theo['omega_y rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']**2), color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], vorticity_theo['omega_z rms']*(settings.Re/wall_UQ_theo['Re_tau_spot']**2), color='tab:green')
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_x rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), color='tab:blue', label_name=r"$\omega_{x~RMS}^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_y rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), color='tab:orange', label_name=r"$\omega_{y~RMS}^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], vorticity['omega_z rms']*(settings.Re/wall_UQ['Re_tau_spot']**2), color='tab:green', label_name=r"$\omega_{z~RMS}^{+}$")
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=2*wall_UQ['Re_tau_spot'])
    fig.chg_y_axis('',axis_low_bound=0,axis_high_bound=0.5)
    fig.custom_layout(enableLegend=True, direction='in')
    fig.add_ticks_x()
    fig.add_ticks_y()
    fig.display()
    
    
def read_and_plot_reynolds_stresses(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the reynolds stresses.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_vel_stats_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/ReynoldsStresses/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/ReynoldsStresses/'
    
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    ####################################################
    #dUfUf/dt
    reynolds_stresses = {
        'name': 'Reynolds_stresses',
            
        'Wall normal coordinate': None,
        'uu': None,
        'vv': None,
        'ww': None,
        'uv': None,
        'uw': None,
        'vw': None
        }
    
    reynolds_stresses = read_stats(reynolds_stresses, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_UQ_theo = read_U_main_quantities(current_path_vel_stats_theo)
    
        # global Mean_RMS
        reynolds_stresses_theo = {
            'name': 'Reynolds_stresses',
                
            'Wall normal coordinate': None,
            'uu': None,
            'vv': None,
            'ww': None,
            'uv': None,
            'uw': None,
            'vw': None
            }
        
        reynolds_stresses_theo = read_stats(reynolds_stresses_theo, current_path_postProcess_theo)
    
    # vorticity RMS
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], -reynolds_stresses['uv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$-\overline{uv}^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], -reynolds_stresses['uv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, linestyle=':', linewidth=2.5, color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], -reynolds_stresses_theo['uv']*(settings.Re/wall_UQ_theo['Re_tau_spot'])**2, color='tab:blue')
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], -reynolds_stresses['uv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:blue', label_name=r"$-\overline{uv}^{+}$")
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=2*wall_UQ['Re_tau_spot'])
    fig.chg_y_axis('Reynolds stresses')
    fig.custom_layout(enableLegend=True)
    fig.add_ticks_x()
    fig.add_ticks_y()
    fig.display()
    
    # fig = CFD_plot('full')
    # fig.add_plot(mesh.Yc, reynolds_stresses['uu']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:blue', label_name=r"$\overline{uu}^{+}$")
    # fig.add_plot(mesh.Yc, reynolds_stresses['vv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:red', label_name=r"$\overline{vv}^{+}$")
    # fig.add_plot(mesh.Yc, reynolds_stresses['ww']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:green', label_name=r"$\overline{ww}^{+}$")
    # fig.chg_x_axis(r'$y^{*}$',axis_low_bound=0,axis_high_bound=1)
    # fig.chg_y_axis('Reynolds stresses',axis_low_bound=0,axis_high_bound=8)
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
    
    # fig = CFD_plot('full')
    # fig.add_plot(mesh.Yc, reynolds_stresses['uu']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:blue', label_name=r"$\overline{uu}^{+}$")
    # fig.add_plot(mesh.Yc, reynolds_stresses['vv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:red', label_name=r"$\overline{vv}^{+}$")
    # fig.add_plot(mesh.Yc, reynolds_stresses['ww']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:green', label_name=r"$\overline{ww}^{+}$")
    # fig.chg_x_axis(r'$y^{*}$',axis_low_bound=2*10**(-3),axis_high_bound=1,axis_scale='log')
    # fig.chg_y_axis('Reynolds stresses',axis_low_bound=0,axis_high_bound=8)
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x(axis_scale='log')
    # fig.add_ticks_y()
    # fig.display()
    
    # fig = CFD_plot('full')
    # fig.add_plot(mesh.Yc, -reynolds_stresses['uv']*(settings.Re/wall_UQ['Re_tau_spot'])**2, color='tab:blue', label_name=r"$-\overline{uv}^{+}$")
    # fig.chg_x_axis(r'$y^{*}$',axis_low_bound=0,axis_high_bound=1)
    # fig.chg_y_axis('Reynolds stresses')
    # fig.custom_layout(enableLegend=True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
    

    
    
def read_and_plot_premultiplied_spectrums(settings:Settings, mesh:CFD_mesh, y_plus:int, current_iteration = None):
    """Plot 2D contours using the CFD_plot class.
    
    Plot the premultiplied spectrums.
    First read .h5 file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (nz,ny/2,nx),
                ...
                'variable_n': array of shape (nz,ny/2,nx)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/Spectrums/'
    
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    # read statistic profiles to ascii file
    if (current_iteration is None):
        fnam = 'Spectrums' + '.h5'
    else:
        fnam = 'Spectrums_' + str(current_iteration*1000) + '.h5'
    
    Spectrums = {
        'name': 'Spectrums',
        
        'uu': None,
        'vv': None,
        'ww': None,
        'uv': None
        }
    
    import h5py
    
    f = h5py.File(current_path_postProcess + '/' + fnam, 'r')
    Spectrums['uu'] = np.array(f['uu_spectrum'])
    Spectrums['vv'] = np.array(f['vv_spectrum'])
    Spectrums['ww'] = np.array(f['ww_spectrum'])
    Spectrums['uv'] = np.array(f['uv_spectrum'])
    Spectrums['lambda_x_plus'] = np.array(f['lambda_x_plus'])
    Spectrums['lambda_z_plus'] = np.array(f['lambda_z_plus'])
    f.close()
    
    j_chosen = find_nearest(mesh.Yc*wall_UQ['Re_tau_spot'], y_plus)

    x_1, x_2 = np.meshgrid(np.transpose(Spectrums['lambda_x_plus']),Spectrums['lambda_z_plus'])
    
    alpha1 = 0.125
    alpha2 = 0.625
    
    extent_x = Spectrums['lambda_x_plus'].shape[0]
    extent_z = Spectrums['lambda_z_plus'].shape[0]
    
    real_y_plus = mesh.Yc[j_chosen]*wall_UQ['Re_tau_spot']
    
    #### UU
    
    sliceToPlot = Spectrums['uu'][:extent_z,j_chosen,:extent_x] * ( wall_UQ['Re_tau_spot'] / settings.Re )**2

    fig = CFD_plot('full')
    fig.colorbar = fig.ax.contourf(x_1,x_2,sliceToPlot,levels=[alpha1*np.max(sliceToPlot),alpha2*np.max(sliceToPlot),np.max(sliceToPlot)],colors=('g', 'r'), extend='max')
    fig.chg_x_axis(r'$\lambda_{x}^{+}$',axis_scale='log')
    fig.chg_y_axis(r'$\lambda_{z}^{+}$',axis_scale='log')
    fig.custom_layout()
    fig.add_title(r'$k_{x}^{+} k_{z}^{+} E_{uu}^{+},~y^{+}=$'+str(round(real_y_plus,2)))
    fig.display()
    
    #### VV
    
    sliceToPlot = Spectrums['vv'][:extent_z,j_chosen,:extent_x] * ( wall_UQ['Re_tau_spot'] / settings.Re )**2

    fig = CFD_plot('full')
    fig.colorbar = fig.ax.contourf(x_1,x_2,sliceToPlot,levels=[alpha1*np.max(sliceToPlot),alpha2*np.max(sliceToPlot),np.max(sliceToPlot)],colors=('g', 'r'), extend='max')
    fig.chg_x_axis(r'$\lambda_{x}^{+}$',axis_scale='log')
    fig.chg_y_axis(r'$\lambda_{z}^{+}$',axis_scale='log')
    fig.custom_layout()
    fig.add_title(r'$k_{x}^{+} k_{z}^{+} E_{vv}^{+},~y^{+}=$'+str(round(real_y_plus,2)))
    fig.display()
    
    #### WW
    
    sliceToPlot = Spectrums['ww'][:extent_z,j_chosen,:extent_x] * ( wall_UQ['Re_tau_spot'] / settings.Re )**2

    fig = CFD_plot('full')
    fig.colorbar = fig.ax.contourf(x_1,x_2,sliceToPlot,levels=[alpha1*np.max(sliceToPlot),alpha2*np.max(sliceToPlot),np.max(sliceToPlot)],colors=('g', 'r'), extend='max')
    fig.chg_x_axis(r'$\lambda_{x}^{+}$',axis_scale='log')
    fig.chg_y_axis(r'$\lambda_{z}^{+}$',axis_scale='log')
    fig.custom_layout()
    fig.add_title(r'$k_{x}^{+} k_{z}^{+} E_{ww}^{+},~y^{+}=$'+str(round(real_y_plus,2)))
    fig.display()
    
    #### UV
    
    sliceToPlot = Spectrums['uv'][:extent_z,j_chosen,:extent_x] * ( wall_UQ['Re_tau_spot'] / settings.Re )**2

    fig = CFD_plot('full')
    fig.colorbar = fig.ax.contourf(x_1,x_2,sliceToPlot,levels=[alpha1*np.max(sliceToPlot),alpha2*np.max(sliceToPlot),np.max(sliceToPlot)],colors=('g', 'r'), extend='max')
    fig.chg_x_axis(r'$\lambda_{x}^{+}$',axis_scale='log')
    fig.chg_y_axis(r'$\lambda_{z}^{+}$',axis_scale='log')
    fig.custom_layout()
    fig.add_title(r'$k_{x}^{+} k_{z}^{+} E_{uv}^{+},~y^{+}=$'+str(round(real_y_plus,2)))
    fig.display()

    
def read_and_plot_premultiplied_spectrums_temperature(settings:Settings, mesh:CFD_mesh, y_plus:int, current_iteration = None):
    """Plot 2D contours using the CFD_plot class.
    
    Plot the premultiplied spectrums.
    First read .h5 file and store values in a dictionnary. Read also main quantities of the flow (Re_tau).
    Then plot the data
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (nz,ny/2,nx),
                ...
                'variable_n': array of shape (nz,ny/2,nx)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
            
    """
    
    # save stats
    current_path_t_stats = settings.path_postProcess+'/TemperatureStats/'
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/TemperatureSpectrums/'
    
    wall_TQ = read_T_main_quantities(current_path_t_stats, current_iteration)
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    # read statistic profiles to ascii file
    if (current_iteration is None):
        fnam = 'Spectrums' + '.h5'
    else:
        fnam = 'Spectrums_' + str(current_iteration*1000) + '.h5'
    
    Spectrums = {
        'name': 'Spectrums',
        
        'tt': None
        }
    
    import h5py
    
    f = h5py.File(current_path_postProcess + '/' + fnam, 'r')
    Spectrums['tt'] = np.array(f['tt_spectrum'])
    Spectrums['lambda_x_plus'] = np.array(f['lambda_x_plus'])
    Spectrums['lambda_z_plus'] = np.array(f['lambda_z_plus'])
    f.close()
    
    j_chosen = find_nearest(mesh.Yc*wall_UQ['Re_tau_spot'], y_plus)

    x_1, x_2 = np.meshgrid(np.transpose(Spectrums['lambda_x_plus']),Spectrums['lambda_z_plus'])
    
    alpha1 = 0.125
    alpha2 = 0.625
    
    extent_x = Spectrums['lambda_x_plus'].shape[0]
    extent_z = Spectrums['lambda_z_plus'].shape[0]
    
    real_y_plus = mesh.Yc[j_chosen]*wall_UQ['Re_tau_spot']
    
    #### TT
    
    sliceToPlot = Spectrums['tt'][:extent_z,j_chosen,:extent_x] * ( 1 / wall_TQ['T_tau_spot'] )**2

    fig = CFD_plot('full')
    fig.colorbar = fig.ax.contourf(x_1,x_2,sliceToPlot,levels=[alpha1*np.max(sliceToPlot),alpha2*np.max(sliceToPlot),np.max(sliceToPlot)],colors=('g', 'r'), extend='max')
    fig.chg_x_axis(r'$\lambda_{x}^{+}$',axis_scale='log')
    fig.chg_y_axis(r'$\lambda_{z}^{+}$',axis_scale='log')
    fig.custom_layout()
    fig.add_title(r'$k_{x}^{+} k_{z}^{+} E_{tt}^{+},~y^{+}=$'+str(round(real_y_plus,2)))
    fig.display()
    
    
def read_and_plot_basic_temperature_stats(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the mean, RMS, skewness and flatness velocities. Plot the diagnostic function, the turbulent Prandtl and the temperature-velocity correlations.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau, T_tau and Nu).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_vel_stats_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/TemperatureStats/'
    
    wall_TQ = read_T_main_quantities(current_path_postProcess, current_iteration)
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    ####################################################
    Mean_Correlations = {
        'name': 'Mean_Correlations',
        
        'Wall normal coordinate': None,
        'T_mean': None,
        'T_rms': None,
        'UfTf': None,
        'VfTf': None,
        'WfTf': None,
        'turb_Pr': None
        }
    
    Mean_Correlations = read_stats(Mean_Correlations, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_TQ_theo = read_T_main_quantities(current_path_postProcess_theo)
        wall_UQ_theo = read_U_main_quantities(current_path_vel_stats_theo)
    
        Mean_Correlations_theo = {
            'name': 'Mean_Correlations',
            
            'Wall normal coordinate': None,
            'T_mean': None,
            'T_rms': None,
            'UfTf': None,
            'VfTf': None,
            'WfTf': None,
            'turb_Pr': None
            }
        
        Mean_Correlations_theo = read_stats(Mean_Correlations_theo, current_path_postProcess_theo)
        
    ######## mean
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_mean']/wall_TQ['T_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_mean']/wall_TQ['T_tau_spot'], linestyle=':', linewidth=2.5, color='k')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_Correlations_theo['T_mean']/wall_TQ_theo['T_tau_spot'])
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_mean']/wall_TQ['T_tau_spot'])
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=5*10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
        fig.chg_y_axis(r'$\overline{\theta}^{+}$', axis_low_bound=0,axis_high_bound=40)
    else:
        fig.chg_y_axis(r'$\overline{T}^{+}$', axis_low_bound=-30,axis_high_bound=0)
    fig.custom_layout()
    fig.add_ticks_x(axis_scale='log')
    fig.add_ticks_y()
    fig.display()

    ######## RMS
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/wall_TQ['T_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/wall_TQ['T_tau_spot'], linestyle=':', linewidth=2.5, color='k')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_Correlations_theo['T_rms']/wall_TQ_theo['T_tau_spot'])
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/wall_TQ['T_tau_spot'])
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
        fig.chg_y_axis(r"$\theta^{'+}_{RMS}$")
    else:
        fig.chg_y_axis(r"$T^{'+}_{RMS}$")
    fig.custom_layout()
    fig.add_ticks_x()
    fig.add_ticks_y()
    fig.display()
    
    # ######## diagnostic function
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
    #         fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/Mean_Correlations['T_mean'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k', label_name=r"$\frac{\overline{\theta^{'}_{RMS}}}{\overline{\theta}}$")
    #     else:
    #         fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/Mean_Correlations['T_mean'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k', label_name=r"$\frac{\overline{T^{'}_{RMS}}}{\overline{T}}$")
    #         fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/Mean_Correlations['T_mean'], linestyle=':', linewidth=2.5, color='k')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_Correlations_theo['T_rms']/Mean_Correlations_theo['T_mean'])
    # else:
    #     if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
    #         fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/Mean_Correlations['T_mean'], label_name=r"$\frac{\overline{\theta^{'+}_{RMS}}}{\overline{\theta}}$")
    #     else:
    #         fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['T_rms']/Mean_Correlations['T_mean'], label_name=r"$\frac{\overline{T^{'+}_{RMS}}}{\overline{T}}$")
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_y_axis('',axis_low_bound=-5,axis_high_bound=5)
    # fig.custom_layout(enableLegend = True)
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()

    ######## temperature_velocity correlations
    fig = CFD_plot('full')
    if (plotWithTheo):
        if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['UfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$\overline{u^{'} \theta^{'}}^{+}$")
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['VfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$\overline{v^{'} \theta^{'}}^{+}$")
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['WfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\overline{w^{'} \theta^{'}}^{+}$")
        else:
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['UfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$\overline{u^{'} T^{'}}^{+}$")
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['VfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$\overline{v^{'} T^{'}}^{+}$")
            fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['WfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:green', label_name=r"$\overline{w^{'} T^{'}}^{+}$")
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['UfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], linestyle=':', linewidth=2.5, color='tab:blue')
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['VfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], linestyle=':', linewidth=2.5, color='tab:orange')
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['WfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], linestyle=':', linewidth=2.5, color='tab:green')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], (Mean_Correlations_theo['UfTf']/wall_TQ_theo['T_tau_spot'])*settings.Re/wall_UQ_theo['Re_tau_spot'], color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], (Mean_Correlations_theo['VfTf']/wall_TQ_theo['T_tau_spot'])*settings.Re/wall_UQ_theo['Re_tau_spot'], color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], (Mean_Correlations_theo['WfTf']/wall_TQ_theo['T_tau_spot'])*settings.Re/wall_UQ_theo['Re_tau_spot'], color='tab:green')
    else:
        if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['UfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:blue', label_name=r"$\overline{u^{'} \theta^{'}}^{+}$")
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['VfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:orange', label_name=r"$\overline{v^{'} \theta^{'}}^{+}$")
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['WfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:green', label_name=r"$\overline{w^{'} \theta^{'}}^{+}$")
        else:
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['UfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:blue', label_name=r"$\overline{u^{'} T^{'}}^{+}$")
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['VfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:orange', label_name=r"$\overline{v^{'} T^{'}}^{+}$")
            fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], (Mean_Correlations['WfTf']/wall_TQ['T_tau_spot'])*settings.Re/wall_UQ['Re_tau_spot'], color='tab:green', label_name=r"$\overline{w^{'} T^{'}}^{+}$")
    fig.chg_x_axis(r'$y^{+}$', axis_low_bound=0, axis_high_bound=wall_UQ['Re_tau_spot'])
    fig.chg_y_axis('')
    fig.custom_layout(enableLegend = True)
    fig.add_ticks_x()
    fig.add_ticks_y()
    fig.display()
                    
    ######## Turb Prandtl
    # fig = CFD_plot('full')
    # fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Mean_Correlations['turb_Pr'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
    # if (plotWithTheo):
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Mean_Correlations_theo['turb_Pr'])
    # fig.chg_x_axis(r'$y^{+}$', axis_low_bound=0, axis_high_bound=wall_UQ['Re_tau_spot'])
    # fig.chg_y_axis(r"$Pr_{t}$", axis_low_bound=0, axis_high_bound=1.2)
    # fig.custom_layout()
    # fig.display()
        
    ####################################################
    Skewness_Flatness = {
        'name': 'Skewness_Flatness',
        
        'Wall normal coordinate': None,
        'T_skewness': None,
        'T_flatness': None
        }
    
    Skewness_Flatness = read_stats(Skewness_Flatness, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        Skewness_Flatness_theo = {
            'name': 'Skewness_Flatness',
            
            'Wall normal coordinate': None,
            'T_skewness': None,
            'T_flatness': None
            }
        
        Skewness_Flatness_theo = read_stats(Skewness_Flatness_theo, current_path_postProcess_theo)
    
    print(current_path_postProcess)
    
    # ######## Skewness
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_skewness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_skewness'], linestyle=':', linewidth=2.5, color='k')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['T_skewness'])
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_skewness'])
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
    #     fig.chg_y_axis(r"$S(\theta^{'+})$")
    # else:
    #     fig.chg_y_axis(r"$S(T^{'+})$")
    # fig.custom_layout()
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
    
    # ######## Flatness
    # fig = CFD_plot('full')
    # if (plotWithTheo):
    #     fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_flatness'], 5, marker_size=10, marker='s', color='None', markeredgecolor='k')
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_flatness'], linestyle=':', linewidth=2.5, color='k')
    #     fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], Skewness_Flatness_theo['T_flatness'])
    # else:
    #     fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], Skewness_Flatness['T_flatness'])
    # fig.chg_x_axis(r'$y^{+}$',axis_low_bound=0,axis_high_bound=wall_UQ['Re_tau_spot'])
    # if ('Ka' in settings.current_path) or ('PSK' in settings.current_path):
    #     fig.chg_y_axis(r"$F(\theta^{'+})$")
    # else:
    #     fig.chg_y_axis(r"$F(T^{'+})$",axis_low_bound=0,axis_high_bound=6)
    # fig.custom_layout()
    # fig.add_ticks_x()
    # fig.add_ticks_y()
    # fig.display()
        
def read_and_plot_temperature_budgets(settings:Settings, mesh:CFD_mesh, current_iteration = None, plotWithTheo = False):
    """Plot 1D statistics using the CFD_plot class.
    
    Plot the temperature budgets.
    First read .dat file and store values in a dictionnary. Read also main quantities of the flow (Re_tau, T_tau and Nu).
    Then plot the data, with the theoretical one (if the boolean is set to True)
    
        The dictionnaries are formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
    
    Parameters
    ----------
    settings : Settings
        Settings of the simulation
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    current_iteration : Integer
        Current iteration of the field to plot
    plotWithTheo : Boolean
        Trigger the option to add the theoretical data on the plot
            
    """
    
    # save stats
    current_path_vel_stats = settings.path_postProcess+'/VelocityStats/'
    current_path_vel_stats_theo = settings.root_postproc + '/' + settings.theo_path + '/VelocityStats/'
    current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
    current_path_postProcess_theo = settings.root_postproc + '/' + settings.theo_path + '/TemperatureStats/'
    
    wall_TQ = read_T_main_quantities(current_path_postProcess, current_iteration)
    wall_UQ = read_U_main_quantities(current_path_vel_stats, current_iteration)
    
    ####################################################
    #dUfUf/dt
    t_budget = {
        'name': 'TemperatureBudgetsSpot',
        'Wall normal coordinate': None,
        'Advection': None,
        'Turbulent Transport': None,
        'Production': None,
        'Molecular diffusion': None,
        'Dissipation': None
        }
    
    t_budget = read_stats(t_budget, current_path_postProcess, current_iteration)
    
    if (plotWithTheo):
    
        wall_TQ_theo = read_T_main_quantities(current_path_postProcess_theo)
        wall_UQ_theo = read_U_main_quantities(current_path_vel_stats_theo)
    
        # global Mean_RMS
        t_budget_theo = {
            'name': 'TemperatureBudgetsSpot',
            'Wall normal coordinate': None,
            'Advection': None,
            'Turbulent Transport': None,
            'Production': None,
            'Molecular diffusion': None,
            'Dissipation': None
            }
        
        t_budget_theo = read_stats(t_budget_theo, current_path_postProcess_theo)
    
    # plot
    fig = CFD_plot('full')
    if (plotWithTheo):
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Production']          * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:blue', label_name=r"$P^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Turbulent Transport'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:orange', label_name=r"$T^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Molecular diffusion'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:red', label_name=r"$D^{+}$")
        fig.add_scatter(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Dissipation']         * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), 5, marker_size=10, marker='s', color='None', markeredgecolor='tab:purple', label_name=r"$\epsilon^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Production']          * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), linestyle=':', linewidth=2.5, color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Turbulent Transport'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), linestyle=':', linewidth=2.5, color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Molecular diffusion'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), linestyle=':', linewidth=2.5, color='tab:red')
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Dissipation']         * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), linestyle=':', linewidth=2.5, color='tab:purple')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], t_budget_theo['Production']          * (settings.Re / (wall_TQ_theo['T_tau_spot']*wall_UQ_theo['Re_tau_spot'])**2), color='tab:blue')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], t_budget_theo['Turbulent Transport'] * (settings.Re / (wall_TQ_theo['T_tau_spot']*wall_UQ_theo['Re_tau_spot'])**2), color='tab:orange')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], t_budget_theo['Molecular diffusion'] * (settings.Re / (wall_TQ_theo['T_tau_spot']*wall_UQ_theo['Re_tau_spot'])**2), color='tab:red')
        fig.add_plot(mesh.Yc*wall_UQ_theo['Re_tau_spot'], t_budget_theo['Dissipation']         * (settings.Re / (wall_TQ_theo['T_tau_spot']*wall_UQ_theo['Re_tau_spot'])**2), color='tab:purple')
    else:
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Production']          * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), color='tab:blue', label_name=r"$P^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Turbulent Transport'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), color='tab:orange', label_name=r"$T^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Molecular diffusion'] * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), color='tab:red', label_name=r"$D^{+}$")
        fig.add_plot(mesh.Yc*wall_UQ['Re_tau_spot'], t_budget['Dissipation']         * (settings.Re / (wall_TQ['T_tau_spot']*wall_UQ['Re_tau_spot'])**2), color='tab:purple', label_name=r"$\epsilon^{+}$")
    fig.chg_x_axis(r'$y^{+}$',axis_low_bound=10**(-1),axis_high_bound=wall_UQ['Re_tau_spot'], axis_scale='log')
    fig.chg_y_axis('')
    fig.custom_layout(enableLegend=True)
    fig.add_title(r'$T$' + ' budget')
    fig.display()