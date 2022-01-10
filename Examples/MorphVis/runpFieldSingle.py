import sys
import time
from math import *
import numpy as np
import os
#first create the seed points of the tessellation
#command = "./build/setCentres 1.0 >output"
#os.system(command)
#now move the seed points file to the logsMorph directory
command = "cp centres.inp ./logsMorph"
os.system(command)
command = "cp pFieldSingle.json ./logsMorph"
os.system(command)
#now run the main program
command = "./build/pFieldSingle pFieldSingle.json "
os.system(command)
print("script finished")
