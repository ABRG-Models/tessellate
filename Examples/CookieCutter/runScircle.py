import sys
import time
from math import *
import numpy as np
import os
command = "cp sCircle.json logsSingleCircle"
os.system("command")
for j in range(0,1):
    command = "./build/setCentres 1.0 >output"
    os.system(command)
    command = "cp centres.inp logsSingleCircle"
    os.system(command)
    command = "cp centres.inp logsSingleCircle/centres" + str(j) + ".inp"
    os.system(command)
    command = "./build/singleCircle sCircle.json " + str(j) + " >output"
    os.system(command)
    time.sleep(10)
print("script finished");
