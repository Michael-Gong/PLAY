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
import random

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


length = 1.0 # 100000*wavelength
ni     = 50*denunit
qi     = q0*79
mi     = m0*1836*189
IS_factor = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)

gg       = np.linspace(1,10010,1001) # Energy_e only for 0.51MeV~5.1GeV
dN_dldw  = np.zeros([1001,5001]) 
dN_dl_sample = np.zeros([1001,5001]) 
for i in range(1001):
    en_photon = np.linspace(1e-4,0.5,5001)*gg[i]*m0*v0**2 # Energy_ph for 0.1~gg/2
    for j in range(5001):
	      dN_dldw[i,j] = np.log(2*gg[i]*(gg[i]*m0*v0**2-en_photon[j])/en_photon[j])/en_photon[j]*(en_photon[-1]-en_photon[-2])
print('dN_dl is constructed !')
dN_dl_sumph = np.sum(dN_dldw,axis=1)
print('dN_dl_sumph is constructed !')
for i in range(1001):
    for j in range(5001):
	      dN_dl_sample[i,j] = np.sum(dN_dldw[i,0:j],axis=0)/dN_dl_sumph[i]
print('dN_dl_sample is constructed !')

np.savetxt('./grid_x_gg.txt',gg)
np.savetxt('./grid_y_ph.txt',np.linspace(1e-4,0.5,5001))
np.savetxt('./dN_dl.txt',dN_dldw*IS_factor)
np.savetxt('./dN_dl_sumph.txt',dN_dl_sumph*IS_factor)
np.savetxt('./dN_dl_sample.txt',dN_dl_sample)
 
table_1 = dN_dl_sumph*IS_factor
table_2 = dN_dl_sample
x_gg    = gg
y_ph    = np.linspace(1e-4,0.5,5001)

e_px = 1000.0
e_gg = (1+e_px**2.0)**0.5
e_vx = e_px/e_gg*v0
e_x  = 0.0
material_l  = 1.0
ni     = 50*denunit
qi     = q0*79
mi     = m0*1836*189
IS_factor = ni*q0**4*qi**2/(4*pi*epsilon0)**3*(16.0/(3*v0**3*h_planck/2/pi))/(m0**2*v0**2)

emission_depth = 0.0
t_time = 0.0
dt     = 3.3e-15 # s
t_stop = 3.3e-12  # s
say_out= 1000.0
step_i = 0
while t_time < t_stop:
    step_i = step_i + 1 
    t_time = t_time + dt
    e_x = e_x + e_vx*dt
    find_i = np.min(np.where(e_gg < table_1))
    rat_1  = (e_gg-x_gg[find_i-1])/(gg[find_i]-gg[find_i-1])
    rat_2  = (x_gg[find_i]-e_gg)/(gg[find_i]-gg[find_i-1])
    emission_depth = emission_depth + e_vx*dt*(table_1[find_i-1]*rat_2+table_1[find_i]*rat_1)
    print('Emission_depth: ',emission_depth)
    if (emission_depth >= 1.0):
        emission_depth = 0.0
        monte_carlo = random.uniform(0, 1)
        find_j = np.min(np.where(monte_carlo < table_2[find_i,:]))
        print('Monte_carlo:',monte_carlo,'; find_j:',find_j)
        rat_3  = (monte_carlo-table_2[find_i,find_j-1])/(table_2[find_i,find_j]-table_2[find_i,find_j-1])
        rat_4  = (table_2[find_i,find_j]-monte_carlo)/(table_2[find_i,find_j]-table_2[find_i,find_j-1])
        print('y_ph[find_j-1]: ',y_ph[find_j-1],'; y_ph[find_j]: ',y_ph[find_j],'; rat_3: ',rat_3,'; rat_4: ',rat_4)
        print('y_ph: ',(y_ph[find_j-1]*rat_4+y_ph[find_j]*rat_3))
        photon_en = (y_ph[find_j-1]*rat_4+y_ph[find_j]*rat_3)*e_gg
        #e_gg = e_gg - photon_en #radiation recoil
        print('Emitting photon with energy: ',photon_en*0.51,' MeV')
    if ( int(step_i)%10 == 9 ):
        print('finished time:',t_time)

