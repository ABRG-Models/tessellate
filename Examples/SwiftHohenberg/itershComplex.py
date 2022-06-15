import sys
import time
from math import *
import numpy as np
import os
command = "cp shComplex.json logsSwiftHohenberg"
os.system("command")
for j in range(0,10):
    command = "./build/shComplex  shComplex.json >output" + str(j)
    os.system(command)
    time.sleep(10)
print("script finished");
