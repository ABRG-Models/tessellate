import sys
import time
from math import *
import numpy as np
import os
command = "cp shComplex.json ./logsTemp"
os.system("command")
for j in range(0,1):
    command = "./build/showPinwheel  shComplex.json >output"
    os.system(command)
    time.sleep(10)
print("script finished");
