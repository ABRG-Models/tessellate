import sys
import time
from math import *
import numpy as np
import os
removeFiles = "rm logsDisplayCircles/*.data logsDisplayCircles/*.txt logsDisplayCircles/*.out logsDisplayCircles/*.h5 logsDisplayCircles/*.png logsDisplayCircles/*.mat"
#os.system(removeFiles)
command = "cp dCircles.json displayCircles.cpp logsDisplayCircles"
os.system(command)
for j in range(0,2):
    command = "./build/setCentres 1.0 >output"
    os.system(command)
    command = "cp centres.inp logsDisplayCircles"
    os.system(command)
    command = "cp centres.inp logsDisplayCircles/centres" + str(j) + ".inp"
    os.system(command)
    command = "sh ./dCirclesjson.sh 36.0 0.0 12.0"
    os.system(command)
    command = "./build/displayCircles dCircles.json " + str(j) + " >output"
    os.system(command)
    time.sleep(10)
print("script finished");
