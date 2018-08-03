#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from scipy.integrate import quad
  
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
  wavelength=     1.00e-6
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

  font_size = 20
  
  to_path='./'
## bremsstrhlung spectrum ##
  gg     = 1000.0
  length = 1.0 # 100000*wavelength
  ni     = 50*denunit
  qi     = q0*79
  mi     = m0*1836*189
  en_photon = np.logspace(-5,np.log10(0.5),1000)*gg*m0*v0**2
  dI_dEdl   = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)*np.log(2*gg*(gg*m0*v0**2-en_photon)/en_photon)*length
  
  plt.subplot(2,2,1)
  plt.plot(en_photon/(1.6e-13),dI_dEdl,'-b',linewidth=4, label=r'$\frac{dW}{dld\hbar\omega}$'+'for $\gamma=1000$',zorder=0)
  #### manifesting colorbar, changing label and axis properties ####
  #plt.xlim(0,520)
  #plt.ylim(5,45)
  plt.xlabel('Photon Energy [MeV]',fontdict=font)
  plt.ylabel(r'$\frac{dW}{dld\hbar\omega} [m^{-1}]$',fontdict=font)
  plt.xscale('log')
  plt.yscale('log')
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=15,framealpha=1.0)

  gg = np.logspace(0,3,1000)
  dI_dl     = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)*gg/2.0*m0*v0**2*np.log(4.0*gg)*length
  averaged_energy = np.log(4*gg)

  plt.subplot(2,2,2)
  plt.plot(gg, dI_dl,'-k',linewidth=4, label=r'$\frac{dW}{dl}$'+'for $E_{ave}=ln(4\gamma)$',zorder=0)
#  plt.plot(gg, dI_dl,'-k',linewidth=4, label=r'$\frac{dW}{dl}$',zorder=0)
  #### manifesting colorbar, changing label and axis properties ####
  #plt.xlim(0,520)
  #plt.ylim(5,45)
  plt.xlabel('$\gamma$',fontdict=font)
  plt.ylabel(r'$\frac{dW}{dl} [J*m^{-1}]$',fontdict=font)
  plt.xscale('log')
  plt.yscale('log')
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=15,framealpha=1.0)



  gg     = 1000.0
  length = 1.0 # 100000*wavelength
  ni     = 50*denunit
  qi     = q0*79
  mi     = m0*1836*189
  en_photon = np.linspace(0.0001,0.5,5001)*gg*m0*v0**2
  dN_dEdl   = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)*np.log(2*gg*(gg*m0*v0**2-en_photon)/en_photon)*length*1.6e-13/en_photon
  plt.subplot(2,2,3)
  plt.plot(en_photon/(1.6e-13),dN_dEdl,'-b',linewidth=4, label=r'$\frac{dN}{dld\hbar\omega}$'+'for $\gamma=1000$',zorder=0)
  #### manifesting colorbar, changing label and axis properties ####
  #plt.xlim(0,520)
  #plt.ylim(5,45)
  plt.xlabel('Photon Energy [MeV]',fontdict=font)
  plt.ylabel(r'$\frac{dN}{dld\hbar\omega} [MeV^{-1}m^{-1}]$',fontdict=font)
  plt.xscale('log')
  plt.yscale('log')
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=15,framealpha=1.0)


  length = 1.0 # 100000*wavelength
  ni     = 50*denunit
  qi     = q0*79
  mi     = m0*1836*189
  gg     = np.logspace(0,3,1001)
  dN_dl  = np.zeros([1001,5001])
  for i in range(1001):
      en_photon = np.linspace(0.0001,0.5,5001)*gg[i]*m0*v0**2
      for j in range(5001):
	      dN_dl[i,j] = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)*np.log(2*gg[i]*(gg[i]*m0*v0**2-en_photon[j])/en_photon[j])*length*1.6e-13/en_photon[j]*(en_photon[-1]-en_photon[-2])/1.6e-13

  plt.subplot(2,2,4)
  plt.plot(gg,np.sum(dN_dl,axis=1),'-k',linewidth=4, label=r'$\frac{dN}{dl}$',zorder=0)
  #### manifesting colorbar, changing label and axis properties ####
  #plt.xlim(0,520)
  #plt.ylim(5,45)
  plt.xlabel('$\gamma$',fontdict=font)
  plt.ylabel(r'$\frac{dN}{dl} [m^{-1}]$',fontdict=font)
  plt.xscale('log')
  plt.yscale('log')
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=15,framealpha=1.0)


  plt.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95,
                wspace=None, hspace=None)
  #plt.text(60,19,r'$\Delta t=\sqrt{\frac{m_ir_0}{2|q|\rho a_0}}$'+' w/o Rel',fontdict=font)


  fig = plt.gcf()
  fig.set_size_inches(20.5, 14)
  fig.savefig('./brem_spectrum_log.png',format='png',dpi=160)
  plt.close("all")
