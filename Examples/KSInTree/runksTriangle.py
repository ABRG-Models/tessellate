import sys
import time
from math import *
import os
# arg1 logpat
# arg2 dt
# arg3 Dn
# arg4 Dchi
# arg5 Dc
# arg6 scale
# arg7 x extent
# arg8 number of time steps
# arg9 printing interval
# arg10 lstart true means start from previous run
# arg11 Lfixedseed if true we are using a fixed seed
# arg12 Lgraphics if true we use graphics
# arg13 LDn if true we adjust the parameters to the area
logpath = "/home/john/Neuroscience/DirichletRD/Threadbeast/sim/KSInTree/logsKSTriangle"
fileRemove = "rm " + logpath + "/*.out " + logpath + "/*.data " + logpath + "/*.txt " + logpath + "/*.png "
os.system(fileRemove)
copyJson = "cp kst.json " + logpath
os.system(copyJson)
copyCpp = "cp ksTriangle.cpp " + logpath
for i in range (0,50):
    json = "./kst.json"
    program = "./build/ksTriangle ./kst.json >output"
    command = "cp " + logpath + "/nnField_00001.png " + logpath + "/nnField1" + str(i) + ".png"
    command1 = "cp " + logpath + "/nnFieldC_00003.png " + logpath + "/nnFieldC" + str(i) + ".png"
    command2 = "cp " + logpath + "/Tesselation1.png " + logpath + "/Tessellation1" + str(i) + ".png"
    print(program)
    print(command)
    print(command1)
    #print(command1)
    os.system(program)
    os.system(command)
    os.system(command1)
    os.system(command2)

