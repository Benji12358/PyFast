"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class CFD_Plot

@description: General class for plotting. This is an instance of matplotlib.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (AutoMinorLocator,LogLocator,NullFormatter,MultipleLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
    	return tuple(i/inch for i in tupl[0])
    else:
       	return tuple(i/inch for i in tupl)
       
params=1
if params==1:   
   
    k=64
    tsize=30
    lsize=tsize+10
    tiksize=tsize+5

    from collections import OrderedDict
    linestyles_dict = OrderedDict(
    [('loosely dotted',        (0, (1, 4))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 3))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (4, 4, 1, 4, 1, 4))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
   
    ls1='solid'#linestyles_dict['solid']
    ls2=linestyles_dict['dashed']
    ls3=linestyles_dict['dotted']
    ls4=linestyles_dict['densely dashed']
    ls5=linestyles_dict['dashdotted']
    ls6=linestyles_dict['densely dashdotted']
    ls7=linestyles_dict['densely dashdotdotted']
    ls8=linestyles_dict['dashdotdotted']
   
    c1='orangered'
    c2='orange'
    c3='forestgreen'
    c4='slategrey'
    c5='royalblue'
    c6='darkred'
    c7='black'
    c8='navy'
   
    m1='o'
    m2='s'
    m3='X'
    m4='D'
    m5='^'
    m6='*'
    m7='v'
    m8='P'
   
    me=5
    ms='10'
    
    plt.rcParams.update({'text.usetex': True, \
                   'axes.labelsize': 40, \
                   'axes.titlesize':40, \
                   'legend.fontsize': 40, \
                   'xtick.labelsize': 40, \
                   'ytick.labelsize': 40})
    
    # params = {'axes.thickness': 0.5}
    # mpl.rcParams.update(params)
    
    # plt.rcParams['font.size'] = 22
    
    # plt.rcParams['xtick.major.pad']='8'
    # plt.rcParams['ytick.major.pad']='8'




class CFD_plot:
    
    def __init__(self, mode, num_plots = None, mydpi = 100, xsize=400.0, ysize=150.0, rows=1, custom_ysize=False):
        
        self.linestyles = linestyles_dict
        
        if (mode=='full'):

            if (custom_ysize):
                self.fig = plt.figure(num=None, figsize=mm2inch(300.0, ysize), dpi=mydpi, linewidth=1)
            else:
                self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 200.0), dpi=mydpi, linewidth=1)
            self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

            for axis in ['top','bottom','left','right']:
                self.ax.spines[axis].set_linewidth(1)
        
        elif (mode=='half_y'):

            self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 100.0), dpi=mydpi, linewidth=3)
            self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

            for axis in ['top','bottom','left','right']:
                self.ax.spines[axis].set_linewidth(1)
        
        elif (mode=='half_x'):

            self.fig = plt.figure(num=None, figsize=mm2inch(150.0, 200.0), dpi=mydpi, linewidth=3)
            self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

            for axis in ['top','bottom','left','right']:
                self.ax.spines[axis].set_linewidth(1)

        elif (mode=='multiple'):

            self.fig = plt.figure(num=None, figsize=mm2inch(xsize, ysize), dpi=mydpi, linewidth=3)
            self.number_plots = num_plots
            self.rows = rows

        elif (mode=='3D'):
    
            self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 200.0), dpi=mydpi, linewidth=1)
            self.ax = self.fig.add_subplot(111, projection="3d")
            
    def update_plot_params(self, size):
        
        params = {'axes.labelsize': size,'axes.titlesize':size, 'legend.fontsize': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
        mpl.rcParams.update(params)


    def create_window(self, index_ax_col, index_ax_row=0):

        self.ax = plt.subplot2grid((self.rows, self.number_plots), (index_ax_row, index_ax_col), rowspan=1, colspan=1)

        for axis in ['top','bottom','left','right']:
            self.ax.spines[axis].set_linewidth(1)

    def switch_window(self, index_axe):

        self.ax = self.fig.axes[index_axe]
        
    def adjust_subplots(self, hspace=0.4, wspace=0.4):
    
        plt.subplots_adjust(wspace=wspace, hspace=hspace)

    def custom_layout(self, enableLegend = False, position = 'center left', direction='out', padding=10):
        
        if (enableLegend):
            if (position=='center right'):

                # Shrink current axis by 20%
                box = self.ax.get_position()
                self.ax.set_position([box.x0 + 0.2*box.width, box.y0, box.width * 0.8, box.height])
                
                self.ax.legend(loc='center right', ncol=1, frameon=False, fancybox=False, facecolor=None, edgecolor=None, framealpha=None, bbox_to_anchor=(-0.05, 0.5))
            else:

            		# Shrink current axis by 20%
                box = self.ax.get_position()
                self.ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                
                self.ax.legend(loc='center left', ncol=1, frameon=False, fancybox=False, facecolor=None, edgecolor=None, framealpha=None, bbox_to_anchor=(1, 0.5))
           
        self.ax.tick_params(axis='x')
        self.ax.tick_params(axis='y')
           
        self.ax.tick_params(which='both',direction=direction, bottom=True, top=False, left=True, right=False)
        self.ax.tick_params(which='both',direction=direction,labelbottom=True, labeltop=False, labelleft=True, labelright=False, pad=padding)
           
        self.ax.tick_params(which='both', width=1)
        self.ax.tick_params(which='major',length=15)
        self.ax.tick_params(which='minor',length=8)

        plt.tight_layout()



    def custom_layout_multiple(self, padding=10, labelsize=25, direction='out', withMinor=True):

        if (direction=='in'):
            
            if (withMinor):
                self.ax.tick_params(which='both', direction='in', width=1, length=5, labelsize=labelsize,
                                    labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                                    bottom=True, top=False, left=True, right=False, pad=18)
               
                self.ax.tick_params(which='major',length=10)
                self.ax.tick_params(which='minor',length=6)
            
        # default ticks direction = out
        else:
            
            self.ax.tick_params(which='both', width=1, length=5, labelsize=labelsize,
                                labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                                bottom=True, top=False, left=True, right=False, pad=padding)
               
            self.ax.tick_params(which='major',length=10)
            self.ax.tick_params(which='minor',length=6)
            
        if not withMinor:
            
            self.ax.tick_params(which='minor',length=0)
        
        
    def add_ticks_x(self, ticks_multiple=None, axis_scale = 'Normal'):
        
        if (axis_scale == 'log'):
            # self.ax.xaxis.set_minor_locator(LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12))
            self.ax.xaxis.set_minor_locator(LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12))
        else:           
        
            if (ticks_multiple is not None):
                self.ax.xaxis.set_minor_locator(MultipleLocator(ticks_multiple))
            else:
                self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        
        
    def add_ticks_y(self, ticks_multiple=None, axis_scale = 'Normal'):
        
        if (axis_scale == 'log'):
            self.ax.yaxis.set_minor_locator(LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),numticks=10))
            # self.ax.yaxis.set_minor_locator(LogLocator())
        elif (axis_scale == 'semlog'):
            self.ax.yaxis.set_minor_locator(LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12))
            self.ax.yaxis.set_minor_locator(LogLocator(base=10.0,subs=(-0.2,-0.4,-0.6,-0.8),numticks=12))
        else:
            if (ticks_multiple is not None):
                self.ax.yaxis.set_minor_locator(MultipleLocator(ticks_multiple))
            else:
                self.ax.yaxis.set_minor_locator(AutoMinorLocator())
                

    def add_plot(self, Xaxis, Yaxis, color='k', linestyle='-', label_name = None, linewidth = 2.5):
    
    	self.ax.plot(Xaxis, Yaxis, color=color, linestyle=linestyle, label=label_name, linewidth=linewidth, rasterized=True)


    def add_scatter(self, Xaxis, Yaxis, step, color='k', markeredgecolor='None', label_name = None, marker = 'o', marker_size = 5, linewidth=2.5, zorder=1):
    
    	self.ax.scatter(Xaxis[::step], Yaxis[::step], c=color, edgecolor=markeredgecolor, label=label_name, marker=marker, s=marker_size**2, linewidth=linewidth, zorder=zorder, rasterized=True)
    
    
    def add_contour(self, ax_1, ax_2, array_2D, cmap, linewidths = 0.5, num_levels = 40, linestyles = '-', dashes=None):
        
        tmp = max(np.max(array_2D), abs(np.min(array_2D)))
        colorbar = np.linspace(-tmp,tmp,num_levels)
        
        if (dashes is not None):
            CS = self.ax.contour(ax_1, ax_2, array_2D, colorbar, cmap=cmap, linewidths=linewidths, linestyles=linestyles, dashes=dashes, rasterized=True)
            
            for c in CS.collections:
                c.set_dashes([(0, dashes)])
        else:
            CS = self.ax.contour(ax_1, ax_2, array_2D, colorbar, cmap=cmap, linewidths=linewidths, linestyles=linestyles, rasterized=True)
    
    
    def add_contourf(self, ax_1, ax_2, array_2D, cmap, num_levels = 40, round_for_label=2):
    
        tmp = max(np.max(array_2D), abs(np.min(array_2D)))
        colorbar = np.linspace(-tmp,tmp,num_levels)
        
        self.m0=round(colorbar[int(0.1*num_levels)],round_for_label)        # colorbar min value
        self.m1=round(colorbar[int(0.3*num_levels)],round_for_label)          # colorbar max value
        self.m2=0.0               # colorbar mid value 1
        self.m3=-self.m1              # colorbar mid value 2
        self.m4=-self.m0               # colorbar mid value 3
        
        self.colorbar = self.ax.contourf(ax_1, ax_2, array_2D, colorbar, cmap=cmap)
        
        for c in self.colorbar.collections:
            c.set_rasterized(True)
    
    
    def chg_x_axis(self, axis_name, axis_low_bound = None, axis_high_bound = None, axis_scale = 'linear', axis_ticks = None, labelpad=5):
    
        self.ax.set_xlabel(axis_name, labelpad=labelpad)

        if (axis_low_bound is not None) and (axis_high_bound is not None):
            self.ax.set_xlim(axis_low_bound, axis_high_bound)
        
        self.ax.set_xscale(axis_scale)

        if (axis_ticks is not None):
            self.ax.set_xticks(axis_ticks)
    
    
    def chg_y_axis(self, axis_name, axis_low_bound = None, axis_high_bound = None, axis_scale = 'linear', axis_ticks = None, labelpad=10):
    
        self.ax.set_ylabel(axis_name, labelpad=labelpad)

        if (axis_low_bound is not None) and (axis_high_bound is not None):
            self.ax.set_ylim(axis_low_bound, axis_high_bound)
        
        self.ax.set_yscale(axis_scale)
    
        if (axis_ticks is not None):
            	self.ax.set_yticks(axis_ticks)
    
    
    def add_title(self, title_name, padding=20):
    
    	plt.title(title_name, pad=padding)
        
    
    def add_colorbar(self, orientation, borderpad=-1.5, title=None, titlepad=None):
        
        if (orientation=='horizontal'):
        
            axins = inset_axes(self.ax,
                        width="100%",  
                        height="5%",
                        loc='upper center',
                        borderpad=borderpad
                       )
            
            for c in self.colorbar.collections:
                c.set_edgecolor("face")
            
            cbar = plt.colorbar(self.colorbar, cax=axins, orientation=orientation)
            cbar.set_ticks([self.m0,self.m1,self.m2,self.m3,self.m4])
            cbar.set_ticklabels([self.m0,self.m1,self.m2,self.m3,self.m4])
            cbar.ax.tick_params(labelsize=25)
            cbar.ax.tick_params(which='major',direction='inout',length=8,width=2)
            if (title is not None) and (titlepad is not None):
                cbar.set_label(title, labelpad=titlepad)
            # cbar.set_ticks([self.m1,self.m2,self.m3])
            # cbar.set_ticklabels([self.m1,self.m2,self.m3])
            
        elif (orientation=='vertical'):
        
            axins = inset_axes(self.ax,
                        width="5%",  
                        height="100%",
                        loc='right',
                        borderpad=borderpad
                       )
            
            for c in self.colorbar.collections:
                c.set_edgecolor("face")
            
            cbar = plt.colorbar(self.colorbar, cax=axins, orientation=orientation)
            cbar.set_ticks([self.m0,self.m1,self.m2,self.m3,self.m4])
            cbar.set_ticklabels([self.m0,self.m1,self.m2,self.m3,self.m4])
            cbar.ax.tick_params(labelsize=25)
            cbar.ax.tick_params(which='major',direction='inout',length=8,width=2)
            if (title is not None) and (titlepad is not None):
                cbar.set_label(title, labelpad=titlepad)
            # cbar.set_ticks([self.m1,self.m2,self.m3])
            # cbar.set_ticklabels([self.m1,self.m2,self.m3])
            
        else:
            print('Orientation not possible, Choose between Vertical and Horizontal')
    
    
    def display(self):
    
        plt.show()
        plt.close()
    
    
    def saveAsPDF(self, path, fig_name):
    
    	plt.savefig(path + fig_name + '.pdf', bbox_inches="tight")


    def saveAsPNG(self, path, fig_name):
    
        plt.savefig(path + fig_name + '.png', bbox_inches="tight")


    def saveAsEPS(self, path, fig_name):
    
        plt.savefig(path + fig_name + '.eps', format='EPS', bbox_inches="tight")