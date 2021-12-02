import sys
from math import *
import numpy as np
import scipy as sp
import scipy.stats as stats
#correlate random ajacent 0
logpath = sys.argv[1]
pC = float(sys.argv[2])
adj0 = np.loadtxt(logpath + "/correlate0.data")
rand0 = np.loadtxt(logpath + "/random_correlate0.data")
#correlate random adjacent 1
adj1 = np.loadtxt(logpath + "/correlate1.data")
rand1 = np.loadtxt(logpath + "/random_correlate1.data")
adj0.sort()
adj1.sort()
rand0.sort()
rand1.sort()
adjSize0 = adj0.size
adjSize1 = adj1.size
randSize0 = rand0.size
randSize1 = rand1.size


countL=0
countM=0
countR=0
for j in range(0,adjSize0):
    if (adj0[j] < -0.9):
        countL = countL + 1
    if ((adj0[j] > -0.5) and (adj0[j] < 0.5)):
        countM = countM + 1
    if (adj0[j] > 0.9):
        countR = countR + 1
countM = float(countM)/ (adjSize0 / (1.0 - pC))
countL = float(countL)/ adjSize0
countR = float(countR)/ adjSize0
print("count L = ", countL, " count M = ", countM, " countR = ", countR)
print("left bimodality noCutAdj = ", float(countL) / countM)
print("right bimodality noCutAdj = ", float(countR) / countM)


countL=0
countM=0
countR=0
for j in range(0,randSize0):
    if (rand0[j] < -0.9):
        countL = countL + 1
    if ((rand0[j] > -0.5) and (rand0[j] < 0.5)):
        countM = countM + 1
    if (rand0[j] > 0.9):
        countR = countR + 1
countM = float(countM)/ (randSize0 / (1.0 - pC))
countL = float(countL)/ randSize0
countR = float(countR)/ randSize0
print("count L = ", countL, " count M = ", countM, " countR = ", countR)
print("left bimodality noCutrand = ", float(countL) / countM)
print("right bimodality noCutrand = ", float(countR) / countM)

countL=0
countM=0
countR=0
for j in range(0,adjSize1):
    if (adj1[j] < -0.9):
        countL = countL + 1
    if ((adj1[j] > -0.5) and (adj1[j] < 0.5)):
        countM = countM + 1
    if (adj1[j] > 0.9):
        countR = countR + 1
countM = float(countM)/ (adjSize1 / (1.0 - pC))
countL = float(countL)/ adjSize1
countR = float(countR)/ adjSize1
print("count L = ", countL, " count M = ", countM, " countR = ", countR)
print("left bimodality CutAdj = ", float(countL) / countM)
print("right bimodality CutAdj = ", float(countR) / countM)

countL=0
countM=0
countR=0
for j in range(0,randSize1):
    if (rand1[j] < -0.9):
        countL = countL + 1
    if ((rand1[j] > -0.5) and (rand1[j] < 0.5)):
        countM = countM + 1
    if (rand1[j] > 0.9):
        countR = countR + 1
countM = float(countM)/ (randSize1 / (1.0 - pC))
countL = float(countL)/ randSize1
countR = float(countR)/ randSize1
print("count L = ", countL, " count M = ", countM, " countR = ", countR)
print("left bimodality Cutrand = ", float(countL) / countM)
print("right bimodality Cutrand = ", float(countR) / countM)
