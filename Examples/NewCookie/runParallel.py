import sys
import time
from math import *
import numpy as np
import os

#make sure setCentres is compiled
#command = "./build/setCentres 1.0 >output "
#os.system(command)
command = "cp ./runParallel.py ./logsdCircles/runParallel.py"
os.system(command)

# Create a list of the commands to be run externally
cmdlist=[]
for j in [2]:
    valTemp = 3.6 - (j)*0.5
    Dn = exp(valTemp)
    Dc = Dn*0.3
    seedInc = 53936398
    for k in [1]:
        seed = seedInc*(85*k+j+1)
        Dchi = (1.0 + 0.143 * k) * Dn
        logpath = "./logsdCircles/Dn" + str(j) + str(k) + "s" 
        command = "mkdir " + logpath
        os.system(command)
        command = "./build/setCentres 1.0"
        os.system(command)
        command = "cp ./centres.inp " + logpath
        os.system(command)
        command = "sh ./dCirclesjson.sh " + str(Dn) + " " + str(Dchi) + " " + str(Dc) + " " + logpath
        os.system(command)
        command = "cp dCircles.json " + logpath
        os.system(command)
        hcommand = "./build/dCircles " + logpath + "/dCircles.json 0  > "  + logpath + "/output"
        print ('Full command: {0}'.format(hcommand))
        cmdlist.append(hcommand)

pl = Pool(4) # how many concurrent commands at a time?
for i, returncode in enumerate(pl.imap(partial(call, shell=True), cmdlist)):
    if returncode != 0:
       print("%d command failed: %d" % (i, returncode))

print("script finished")
