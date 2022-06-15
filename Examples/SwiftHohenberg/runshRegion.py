import sys
import time
from math import *
import numpy as np
import os
command = "cp shComplex.json logsSwiftHohenberg"
#first create the seed points of the tessellation
command = "./build/setCentres 1.0 >output"
os.system(command)
#now move the seed points file to the logsMorph directory
command = "cp centres.inp ./logsSwiftHohenberg"
os.system(command)
#now run the main program
os.system("command")
for j in range(0,1):
    command = "./build/shRegion  shComplex.json >output"
    os.system(command)
    time.sleep(10)
print("script finished");
