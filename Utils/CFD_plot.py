"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class CFD_Plot

@description: General class for plotting. This is an instance of matplotlib.
"""

def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
    	return tuple(i/inch for i in tupl[0])
    else:
       	return tuple(i/inch for i in tupl)
       
params=1
if params==1:   
   
    k=64
    tsize=20
    lsize=tsize+7
    tiksize=tsize+5

    from collections import OrderedDict
    linestyles_dict = OrderedDict(
    [('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
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



import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator)
import numpy as np

class CFD_plot:
    
    def __init__(self, mode, num_plots = None, mydpi = 100):
        
        if (mode=='full'):

            self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 200.0), dpi=mydpi, linewidth=3)
            self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

            for axis in ['top','bottom','left','right']:
                self.ax.spines[axis].set_linewidth(2)
        
        if (mode=='half'):

            self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 100.0), dpi=mydpi, linewidth=3)
            self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

            for axis in ['top','bottom','left','right']:
                self.ax.spines[axis].set_linewidth(2)

        elif (mode=='multiple'):

            self.fig = plt.figure(num=None, figsize=mm2inch(300.0, 100.0), dpi=mydpi, linewidth=3)
            self.number_plots = num_plots


    def create_window(self, index_axe):

        self.ax = plt.subplot2grid((1, self.number_plots), (0, index_axe), rowspan=1, colspan=1)

        for axis in ['top','bottom','left','right']:
            self.ax.spines[axis].set_linewidth(2)

    def switch_window(self, index_axe):

        self.ax = self.fig.axes[index_axe]

    def custom_layout(self, enableLegend = False):

		# Shrink current axis by 20%
        box = self.ax.get_position()
        self.ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        if (enableLegend):
            	self.ax.legend(loc='center left',fontsize=tsize, ncol=1, frameon=False, fancybox=False, facecolor=None, edgecolor=None, framealpha=None, bbox_to_anchor=(1, 0.5))
           
        self.ax.tick_params(axis='x', labelsize=tiksize)
        self.ax.tick_params(axis='y', labelsize=tiksize)
           
        self.ax.tick_params(which='both',bottom=True, top=True, left=True, right=True)
        self.ax.tick_params(which='both',labelbottom=True, labeltop=False, labelleft=True, labelright=False)
           
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.tick_params(which='both', width=2)
        self.ax.tick_params(which='major',length=15)
        self.ax.tick_params(which='minor',length=8)

        plt.tight_layout()



    def custom_layout_multiple(self, x_limits, y_limits):

        for index_axe in range(self.number_plots): 
            self.fig.axes[index_axe].tick_params(axis='x', labelsize=tiksize)
            self.fig.axes[index_axe].tick_params(axis='y', labelsize=tiksize)
            
            if i == 0:
                self.fig.axes[index_axe].tick_params(which='both',bottom=True, top=True, left=True, right=False)
                self.fig.axes[index_axe].tick_params(which='both',labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            
                self.fig.axes[index_axe].set_xlim([x_limits[index_axe][0],x_limits[index_axe][1]+2])
                self.fig.axes[index_axe].set_ylim([y_limits[index_axe][0],y_limits[index_axe][1]])
                
            elif i == len(time_values)-1:
                self.fig.axes[index_axe].tick_params(which='both',bottom=True, top=True, left=False, right=True)
                self.fig.axes[index_axe].tick_params(which='both',labelbottom=True, labeltop=False, labelleft=False, labelright=False)
            
                self.fig.axes[index_axe].set_xlim([x_limits[index_axe][0]-2,x_limits[index_axe][1]])
                self.fig.axes[index_axe].set_ylim([y_limits[index_axe][0],y_limits[index_axe][1]])
                
            else:
                self.fig.axes[index_axe].tick_params(which='both',bottom=True, top=True, left=False, right=False)
                self.fig.axes[index_axe].tick_params(which='both',labelbottom=True, labeltop=False, labelleft=False, labelright=False)
            
                self.fig.axes[index_axe].set_xlim([x_limits[index_axe][0]-2,x_limits[index_axe][1]+2])
                self.fig.axes[index_axe].set_ylim([y_limits[index_axe][0],y_limits[index_axe][1]])
                
                self.fig.axes[index_axe].set_yticks([])

                
            # Shrink current axis by 20%
            box = self.fig.axes[index_axe].get_position()
            self.fig.axes[index_axe].set_position([box.x0, box.y0, box.width * 0.8, box.height])
           
            self.fig.axes[index_axe].xaxis.set_minor_locator(AutoMinorLocator(5))
            self.fig.axes[index_axe].yaxis.set_minor_locator(AutoMinorLocator(5))
            self.fig.axes[index_axe].tick_params(which='both', width=2)
            self.fig.axes[index_axe].tick_params(which='major',length=15)
            self.fig.axes[index_axe].tick_params(which='minor',length=8)

        plt.tight_layout()
        plt.subplots_adjust(wspace=0,hspace=0)












        # Shrink current axis by 20%
        box = self.ax.get_position()
        self.ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        if (enableLegend):
                self.ax.legend(loc='center left',fontsize=tsize, ncol=1, frameon=False, fancybox=False, facecolor=None, edgecolor=None, framealpha=None, bbox_to_anchor=(1, 0.5))
           
        self.ax.tick_params(axis='x', labelsize=tiksize)
        self.ax.tick_params(axis='y', labelsize=tiksize)
           
        self.ax.tick_params(which='both',bottom=True, top=True, left=True, right=True)
        self.ax.tick_params(which='both',labelbottom=True, labeltop=False, labelleft=True, labelright=False)
           
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.tick_params(which='both', width=2)
        self.ax.tick_params(which='major',length=15)
        self.ax.tick_params(which='minor',length=8)

        plt.tight_layout()


    def add_plot(self, Xaxis, Yaxis, color='k', linestyle='-', label_name = None):
    
    	self.ax.plot(Xaxis, Yaxis, color=color, linestyle=linestyle, label=label_name)


    def add_scatter(self, Xaxis, Yaxis, step, color='k', markeredgecolor='None', label_name = None, marker = 'o', marker_size = 5):
    
    	self.ax.scatter(Xaxis[::step], Yaxis[::step], c=color, edgecolor=markeredgecolor, label=label_name, marker=marker, s=marker_size**2)
    
    
    def add_contour(self, ax_1, ax_2, array_2D, cmap, linewidths = 0.5):
    
    	tmp = max(np.max(array_2D), abs(np.min(array_2D)))
    	colorbar = np.linspace(-tmp,tmp,100)
    
    	self.ax.contour(ax_1, ax_2, array_2D, colorbar, cmap=cmap, linewidths=linewidths)
    
    
    def add_contourf(self, ax_1, ax_2, array_2D, cmap, num_levels = 40):
    
    	tmp = max(np.max(array_2D), abs(np.min(array_2D)))
    	colorbar = np.linspace(-tmp,tmp,num_levels)
    
    	self.colorbar = self.ax.contourf(ax_1, ax_2, array_2D, colorbar, cmap=cmap)
    
    
    def chg_x_axis(self, axis_name, axis_low_bound = None, axis_high_bound = None, axis_scale = 'linear', axis_ticks = None):
    
        self.ax.set_xlabel(axis_name,fontsize=lsize)

        if (axis_low_bound is not None) and (axis_high_bound is not None):
            self.ax.set_xlim(axis_low_bound, axis_high_bound)
        
        self.ax.set_xscale(axis_scale)

        if (axis_ticks is not None):
            self.ax.set_xticks(axis_ticks)
    
    
    def chg_y_axis(self, axis_name, axis_low_bound = None, axis_high_bound = None, axis_scale = 'linear', axis_ticks = None):
    
        self.ax.set_ylabel(axis_name,fontsize=lsize)

        if (axis_low_bound is not None) and (axis_high_bound is not None):
            self.ax.set_ylim(axis_low_bound, axis_high_bound)
        
        self.ax.set_yscale(axis_scale)
    
        if (axis_ticks is not None):
            	self.ax.set_yticks(axis_ticks)
    
    
    def add_title(self, title_name):
    
    	plt.title(title_name)
        
    
    def add_colorbar(self, orientation):
        
        plt.colorbar(self.colorbar, orientation=orientation)
    
    
    def display(self):
    
        plt.show()
        plt.close()
    
    
    def saveAsPDF(self, path, fig_name):
    
    	plt.savefig(path + fig_name + '.pdf', bbox_inches="tight")


    def saveAsPNG(self, path, fig_name):
    
        plt.savefig(path + fig_name + '.png', bbox_inches="tight")


    def saveAsEPS(self, path, fig_name):
    
        plt.savefig(path + fig_name + '.eps', bbox_inches="tight")