from scipy.integrate import odeint
%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
rc('font',**{'family':'sans-serif','sans-serif':['Palatino']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

font = {'family' : 'helvetica',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
        }

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

x  = np.array([300,700, 1000, 1400, 1800, 2100, 2500, 2800])
y1 = np.array([0.347072628616205,0.411508800940784,0.465992733403173,0.532932605246774,0.562312456714083,0.577617352017823,0.585695805213658,0.589221080550980])
y2 = np.array([0.365215312035121,0.433083091745553,0.490009172579773,0.543862748207501,0.572795375241772,0.585903725526069,0.590774627501026,0.593817542223461])
y3 = np.array([0.392208938352955,0.455703502161627,0.502773937765499,0.564866593650685,0.585020652189350,0.602209045051294,0.6029328818335834,0.605769990061125])


N = 8

ind = np.arange(N) 
width = 100     
#plt.bar(ind, y1, width*0.5, label=r'$User\ density\ =\ 93.3/km^2$', color='b', alpha=0.5, hatch="o")
#plt.bar(ind + 0.5*width, y2, width*0.5, label=r'$User\ density\ =\ 178.3/km^2$', color='g', alpha=0.5, hatch='/') 
#plt.bar(ind+ 1*width, y3, width*0.5, label=r'$User\ density\ =\ 263.1/km^2$', color='r', alpha=0.5, hatch='*')

plt.bar(x- 0.5*width, y1,  width*0.5, label=r'$User\ density\ =\ 93.3/km^2$', color='b', alpha=0.5, hatch="o")
plt.bar(x, y2,  width*0.5, label=r'$User\ density\ =\ 178.3/km^2$', color='g', alpha=0.5, hatch='/') 
plt.bar(x+ 0.5*width, y3,  width*0.5, label=r'$User\ density\ =\ 263.1/km^2$', color='r', alpha=0.5, hatch='*')
        
plt.ylabel(r'$Propotion\ of\ users\ accessing\ the\ LEO-based\ small\ cells$', fontsize=20)
plt.xlabel(r'$Data\ amount\ generated\ by\ each\ user\ (bytes/s)$', fontsize=20)
#plt.xlim(0, 1.0)
plt.ylim(0,0.7)

#plt.title('Scores by group and gender')

#plt.xticks(ind + width / 2, (r'300', r'700', r'1000', r'1400', r'1800', r'2100', r'2500', r'2800'))
plt.legend(loc='best', framealpha=0.0, fontsize=18)
plt.xticks([300,700, 1000, 1400, 1800, 2100, 2500, 2800])

#fig = plt.figure()
#ax2 = fig.add_subplot(111)
#bars = ax2.bar(range(1, 5), range(1, 5), color='yellow', ecolor='black') + \
#    ax2.bar(range(1, 5), [6] * 4, bottom=range(1, 5), color='green', ecolor='black')
#ax2.set_xticks([1.5, 2.5, 3.5, 4.5])

#patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')
#for bar, pattern in zip(bars, patterns):
#    bar.set_hatch(pattern)

plt.tight_layout()
#plt.show()
#plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.02, wspace=0.50)
fig = plt.gcf()
fig.set_size_inches(11, 8)
#fig.set_size_inches(5, 4.5)
fig.savefig('./dibo.png',format='png',dpi=480)
#plt.close("all")
