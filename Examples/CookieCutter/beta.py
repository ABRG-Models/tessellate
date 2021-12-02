import numpy as np
#from sebcolour import Colour as C
#aclr = C.violetred   # Colour for 'adjacent' data
#rclr = C.dodgerblue  # Colour for 'random' data

# Set plotting font defaults
import matplotlib
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)

# Important for svg output of text as 'things that can be edited in inkscape'
import pylab as pl
pl.rcParams['svg.fonttype'] = 'none'

import scipy.special as sc # for the beta *function* used to compute beta dist.

# Compute the beta *distribution*. x is a set of real values in the open range
# (0,1) as an array.
def beta_distribution (x, a, b):
    betadist = (1./sc.beta (a, b)) * np.power(x, a-1) * np.power(1-x, b-1)
    return betadist

x = np.linspace(0.00001, 0.99999, 1000)

fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot (1,1,1);

for alpheqbet in [0.1,0.5,0.99999]:
    bd = beta_distribution (x, alpheqbet, alpheqbet)
    ax.plot (x, bd, label='a=b={0:01f}'.format(alpheqbet))

ax.legend()

pl.show()
