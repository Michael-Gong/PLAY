#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
#import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import multiprocessing as mp
  
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
  wavelength=     0.8e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2 #J_a/lambda^2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 25,  
          }  
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 16,  
          } 
  space_1 = 1
  space_2 = 1
  font_size = 25
  marker_size=0.05
##below is for generating mid transparent colorbar

if __name__ == '__main__':
    to_path = './'

    plt.subplot(1,1,1)
    w0 = 2.4
    y  = np.linspace(-5,5,200)
    Ey = np.exp(-y**2/w0**2)
    plt.plot(y,Ey)
    w0 = w0/1.5
    Ey = np.exp(-(y-1.5)**2/w0**2)+0.5*np.exp(-(y+1.5)**2/w0**2)
    plt.plot(y,Ey,label='d=1.5')
    Ey = np.exp(-(y-2)**2/w0**2)+0.5*np.exp(-(y+2)**2/w0**2)
    plt.plot(y,Ey,label='d=2.0')
    Ey = np.exp(-(y-2.5)**2/w0**2)+0.5*np.exp(-(y+2.5)**2/w0**2)
    plt.plot(y,Ey,label='d=2.5')

        
#    plt.xlim(-7,7)    
#    plt.xlim(-2,28)
    plt.xlabel('$y$ [$\mu m$]',fontdict=font)
    plt.ylabel('$E_y$',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
#    plt.text(20.,5.5,'n='+str(n)+' t='+str(round(time/1.0e-15-7.5*10./3.,2))+' fs',fontdict=font2,color='k')
    plt.legend(fontsize=font_size-5)

    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
    fig = plt.gcf()
    fig.set_size_inches(10,8.5)
    fig.savefig(to_path+'Ey_noGaussian.png',format='png',dpi=160)
    plt.close("all")
    print(to_path+'Ey_noGaussian.png')



