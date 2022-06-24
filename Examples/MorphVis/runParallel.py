import sys
import time
from math import *
import numpy as np
import os

#make sure setCentres is compiled
#command = "./build/setCentres 1.0 >output "
#os.system(command)
command = "cp ./runParallel.py ./logsMorph/runParallel.py"
os.system(command)

# Create a list of the commands to be run externally
cmdlist=[]
for j in range(0,4):
    valTemp = 3.6 - (j)*0.5
    Dn = exp(valTemp)
    Dc = Dn*0.3
    seedInc = 7306026
    for k in range(0,8):
        seed = seedInc*(67*k+j+1)
        Dchi = (1.0 + 0.143 * k) * Dn
        logpath = "./logsMorph/Dn" + str(j) + str(k)
        command = "mkdir " + logpath
        os.system(command)
#        command = "./build/setCentres 1.0"
#        os.system(command)
        command = "cp ./ratCentres.inp " + logpath
        os.system(command)
        command = "sh ./pFieldVisjson.sh " + str(Dn) + " " + str(Dchi) + " " + str(Dc)
        os.system(command)
        command = "cp pFieldVis.json " + logpath
        os.system(command)
        hcommand = "./build/pFieldRem " + logpath + "/pFieldVis.json  0 10000000 " + logpath + " > "  + logpath + "/output"
        print ('Full command: {0}'.format(hcommand))
        cmdlist.append(hcommand)

from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call

pl = Pool(64) # how many concurrent commands at a time?
for i, returncode in enumerate(pl.imap(partial(call, shell=True), cmdlist)):
    if returncode != 0:
       print("%d command failed: %d" % (i, returncode))

print("script finished")
