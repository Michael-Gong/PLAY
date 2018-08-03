#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
 

font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
font_size=20

data = np.loadtxt('./photon_1000.0.txt')
 
# example data
#mu = 100 # mean of distribution
#sigma = 15 # standard deviation of distribution
#x = mu + sigma * np.random.randn(10000)
 
num_bins = 250
# the histogram of the data
n, bins, patches = plt.hist(data, num_bins, facecolor='blue', alpha=0.75)
 
# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
plt.ylabel('Number [MeV$^{-1}$*m$^{-1}$]',fontdict=font)
plt.xlabel('Photon Energy [MeV]',fontdict=font)
plt.title('Spectrum for benchmark')
plt.yscale('log')
 
# Tweak spacing to prevent clipping of ylabel
plt.xticks(fontsize=font_size);
plt.yticks(fontsize=font_size);
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
#plt.legend(loc='best',fontsize=15,framealpha=1.0)


plt.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95,
                wspace=None, hspace=None)
  #plt.text(60,19,r'$\Delta t=\sqrt{\frac{m_ir_0}{2|q|\rho a_0}}$'+' w/o Rel',fontdict=font)
fig = plt.gcf()
fig.set_size_inches(8, 6.5)
fig.savefig('./hist_spectrum.png',format='png',dpi=160)
plt.close("all")
