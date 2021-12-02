import sys
from math import *
import numpy as np
import scipy as sp
import scipy.stats as stats
logpath = sys.argv[1]
#correlate random ajacent 0
adjacent = np.loadtxt(logpath + "/correlate0.data")
random = np.loadtxt(logpath + "/random_correlate0.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.chisquare(histRand, histAdj)
print ("chisquare rand0/corr0 = ", KL)
KL = stats.chisquare(histAdj, histRand)
print ("chisquare corr0/rand0 = ", KL)
#correlate random adjacent 1
adjacent = np.loadtxt(logpath + "/correlate1.data")
random = np.loadtxt(logpath + "/random_correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.chisquare(histRand, histAdj)
print ("chisquare rand1/corr1 = ", KL)
KL = stats.chisquare(histAdj, histRand)
print ("chisquare corr1/rand1 = ", KL)
#correlate adj adj noCut cut
adjacent = np.loadtxt(logpath + "/correlate0.data")
random = np.loadtxt(logpath + "/correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.chisquare(histRand, histAdj)
print ("chisquare corr0/corr1 = ", KL)
KL = stats.chisquare(histAdj, histRand)
print ("chisquare corr1/corr0 = ", KL)
#corr3
adjacent = np.loadtxt(logpath + "/random_correlate0.data")
random = np.loadtxt(logpath + "/random_correlate1.data")
[histAdj,binsAdj] = np.histogram(adjacent,20)
[histRand,binsRand] = np.histogram(random, 20)
KL = stats.chisquare(histRand, histAdj)
print ("chisquare rand0/rand1 = ", KL)
KL = stats.chisquare(histAdj, histRand)
print ("chisquare rand1/rand0 = ", KL)
