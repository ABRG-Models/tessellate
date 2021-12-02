import sys
import time
from math import *
import numpy as np
import os
command = "cp sH.json logsSwiftHohenberg"
os.system("command")
for j in range(0,1):
    command = "./build/setCentres 1.0 >output"
    os.system(command)
    command = "cp centres.inp logsSwiftHohenberg"
    os.system(command)
    command = "cp centres.inp logsSwiftHohenberg/centres" + str(j) + ".inp"
    os.system(command)
    command = "./build/swiftHohenberg  sH.json >output"
    os.system(command)
    time.sleep(10)
print("script finished");
