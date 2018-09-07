#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec


if __name__ == "__main__":
    print ('This is main of module "test2d.py"')
    ######## Constant defined here ########
    pi        =     3.1415926535897932384626
    q0        =     1.602176565e-19 # C
    m0        =     9.10938291e-31  # kg
    v0        =     2.99792458e8    # m/s^2
    kb        =     1.3806488e-23   # J/K
    mu0       =     4.0e-7*np.pi       # N/A^2
    epsilon0  =     8.8541878176203899e-12 # F/m
    h_planck  =     6.62606957e-34  # J s
    wavelength=     1.0e-6
    frequency =     v0*2*pi/wavelength

    exunit    =     m0*v0*frequency/q0
    bxunit    =     m0*frequency/q0
    denunit    =     frequency**2*epsilon0*m0/q0**2
    jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
    print('electric field unit: '+str(exunit))
    print('magnetic field unit: '+str(bxunit))
    print('density unit nc: '+str(denunit))

    font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
    font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 15,  
          }  
    
##below is for generating mid transparent colorbar
    c_red = matplotlib.colors.colorConverter.to_rgba('red')
    c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
    c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
    cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
    cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128) 
##end for transparent colorbar##

##below is for norm colorbar
    class MidpointNormalize(colors.Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            colors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
            # I'm ignoring masked values and all kinds of edge cases to make a
            # simple example...
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####
    gg0 = 1.0
    density_n0 = 1.0*denunit/gg0

    x1   = np.linspace(-10,10,801)*wavelength*gg0
    density_1 = np.zeros([801])
    grid_x1 = x1
    density_1[abs(grid_x1/gg0)<4.0*wavelength] = density_n0

    point_x = 4.0*wavelength
    ex = q0*density_1*(point_x-grid_x1)/(2*epsilon0)/((point_x-grid_x1)+1.e-80)
    ex = sum(ex)*(x1[-1]-x1[-2])
    print('Static frame, Ex:',ex/exunit,'; No Ey')
    

    

    gg0 = 100.0
    density_n0 = 1.0*denunit/gg0

    x1   = np.linspace(-10,10,801)*wavelength*gg0
    density_1 = np.zeros([801])
    grid_x1 = x1
    density_1[abs(grid_x1/gg0)<4.0*wavelength] = density_n0

    point_x = 4.0*wavelength
    ex = q0*density_1*(point_x-grid_x1)/(2*epsilon0)/((point_x-grid_x1)+1.e-80)
    ex = sum(ex)*(x1[-1]-x1[-2])
    print('Lab frame, Ex:',ex/exunit,'; No Ey')
    

    
    
    
