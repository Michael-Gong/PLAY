import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import scipy.integrate as integrate
import scipy.special as special 
from scipy.special import kv

######## Constant defined here ########
pi        =     3.1415926535897932384626
pi2d      =     180./pi
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  

###############read data into array################
part_number=1
nsteps=3999

gg=np.loadtxt('./Data/px_0.txt')
b0=np.loadtxt('./Data/bz_part_0.txt')
t0=np.loadtxt('./Data/t_0.txt')
gg = gg[0]
b0 = b0[1]
t0 = t0[-1]
r0=gg/b0
omega_critic = 1.5*gg**3/r0

circle_number=t0/(2*np.pi*gg/b0)

insert1='./Data/'
insert_n='_0'
photon=np.loadtxt(insert1+'qed_photon'+insert_n+'.txt')
weight_photon=photon
photon = photon*m0*v0**2/(h_planck/6.28*frequency)

def pxpy_to_energy(gamma, weight):
      binsize = 600
      en_grid = np.linspace(10*omega_critic*(1.0/binsize*0.5),10*omega_critic*(1.0-1.0/binsize*0.5),binsize)
      en_bin  = np.linspace(0,10*omega_critic,binsize+1)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(10*omega_critic/binsize)
        print('number:',np.size(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1])]),'; value:',en_value[i])
      return (en_grid, en_value)

value_x,value_y = pxpy_to_energy(photon,weight_photon)
print(value_y)

print('total circle: ',circle_number)

value_y = value_y/6.28/circle_number

#theory_line = norm_fac*3*(2*grid_omega_x*gg/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.6667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.33333,xi)**2)

#norm_x = matplotlib.colors.Normalize()
plt.subplot(1,2,1)
#plt.plot(grid_omega_x,theory_line,'-r',linewidth=3,label='theoretical equation')
plt.plot(value_x, value_y, '-b',linewidth=3,label='QED_photon')
#cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
#cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
#cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
#plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
#### manifesting colorbar, changing label and axis properties ####
#plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,10*omega_critic)
#plt.ylim(0,8.5e-7)
plt.legend(loc='best',fontsize=16,framealpha=1.0)
#plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
#plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)

plt.subplot(1,2,2)
#plt.plot(grid_omega_x, data_omega, '-k',linewidth=1)
#plt.plot(grid_omega_x,theory_line,'-r',linewidth=3)
#plt.plot(x_line, y_line, ':b',linewidth=3)
#cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
#cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
#cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
#plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
#### manifesting colorbar, changing label and axis properties ####
#plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.plot(value_x, value_y, '-b',linewidth=3,label='QED_photon')
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.xscale('log')
plt.yscale('log')
plt.xlim(0,10*omega_critic)
#plt.ylim(1e-9,10e-7)
#plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
#plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
#plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)

fig = plt.gcf()
fig.set_size_inches(18.0, 8.5)
fig.savefig('./Data/spectral_qed_omega.png',format='png',dpi=160)
plt.close("all")

