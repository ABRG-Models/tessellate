import sys
import time
from math import *
import numpy as np
import os
outfilename1 = './lowAlpha.data'
print(outfilename1)
outfile1 = open(outfilename1,'w')
for j in range(0,32):
    logpath = "./RandomLogs/LowAlpha/lowAlpha" + str(j)
    infilename = logpath + '/correlate1.data'
    infile = open(infilename,'r')
    print(infilename)
    lines = infile.readlines()
    line_no = 0
    for line in lines:
        line_no += 1
        outfile1.write('%s' % line)
        print(line)
    infile.close()
outfile1.close()
print("script finished")
