import sys
from math import *
import numpy as np
import scipy as sp
import scipy.stats as stats
#correlate random ajacent 0
adjacent = np.loadtxt("./correlate0.data")
random = np.loadtxt("./random_correlate0.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.entropy(histRand, histAdj)
print ("KL divergence rand0/corr0 = ", KL)
KL = stats.entropy(histAdj, histRand)
print ("KL divergence corr0/rand0 = ", KL)
#correlate random adjacent 1
adjacent = np.loadtxt("./correlate1.data")
random = np.loadtxt("./random_correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.entropy(histRand, histAdj)
print ("KL divergence rand1/corr1 = ", KL)
KL = stats.entropy(histAdj, histRand)
print ("KL divergence corr1/rand1 = ", KL)
#correlate adj adj noCut cut
adjacent = np.loadtxt("./correlate0.data")
random = np.loadtxt("./correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.entropy(histRand, histAdj)
print ("KL divergence corr0/corr1 = ", KL)
KL = stats.entropy(histAdj, histRand)
print ("KL divergence corr1/corr0 = ", KL)
#corr3
adjacent = np.loadtxt("./random_correlate0.data")
random = np.loadtxt("./random_correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.entropy(histRand, histAdj)
print ("KL divergence rand0/rand1 = ", KL)
KL = stats.entropy(histAdj, histRand)
print ("KL divergence rand1/rand0 = ", KL)
