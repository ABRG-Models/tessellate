import sys
import time
from math import *
import numpy as np
import os
#removeFiles = "rm logsdCircles/*.data logsdCircles/*.txt logsdCircles/*.out logsdCircles/*.h5 logsdCircles/*.png logsdCircles/*.mat"
#os.system(removeFiles)
command = "cp dCircles.json dCircles.cpp logsdCircles"
os.system(command)
for j in range(0,1):
    command = "./build/setCentres 1.0 >output"
    os.system(command)
    command = "cp centres.inp logsdCircles"
    os.system(command)
    command = "cp centres.inp logsdCircles/centres" + str(j) + ".inp"
    os.system(command)
    command = "./build/dCircles dCircles.json " + str(j) + " >output"
    os.system(command)
    command = "cp ./logsdCircles/first.h5 ./logsdCircles/first" + str(j) + ".h5"
    os.system(command)
    command = "cp ./logsdCircles/second.h5 ./logsdCircles/second" + str(j) + ".h5"
    os.system(command)
    command = "cp ./logsdCircles/nnField0.png ./logsdCircles/nnField0" + str(j) + ".png"
    os.system(command)
    command = "cp ./logsdCircles/nnField20.png ./logsdCircles/nnField20" + str(j) + ".png"
    os.system(command)
    command = "cp ./logsdCircles/nnLapl20.png ./logsdCircles/nnLapl20" + str(j) + ".png"
    os.system(command)
    time.sleep(10)
print("script finished");
