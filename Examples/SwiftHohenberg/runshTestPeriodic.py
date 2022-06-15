import sys
import time
from math import *
import numpy as np
import os
command = "cp shTestPeriodic.json logsSwiftHohenberg"
os.system("command")
for j in range(0,1):
    command = "./build/shTestPeriodic  shTestPeriodic.json >output1"
    os.system(command)
    time.sleep(10)
print("script finished");
