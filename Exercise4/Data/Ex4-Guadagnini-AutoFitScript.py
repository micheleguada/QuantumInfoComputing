## Exercise 4 - Quantum Information and Computing 2020/21 ##
# This script executes the fit automatically for all the datafiles
# in the folder by using gnuplot.
#
# AUTHOR: Michele Guadagnini - ID 1230663

import time
import glob
from os import system

gnufile = "AutoFit.gnu"

## begin ##
start_time = time.time()
print("Automated fits")

# list of datafiles in the current folder
txtfiles = []
for dat in glob.glob("*.txt"):
    txtfiles.append(dat)

# resetting the fit.log file
print("  resetting the fit.log file...", end='')
open("fit.log","w").close()
print("  Done.")

## gnuplot part ##
for data in txtfiles:
    system("gnuplot -c "+gnufile+" "+data)
    print("  Fit of data in file "+data+" completed.")
    system("mv fit.log "+data[:-4]+"_fit.log") 
            
print("Time elapsed from start [s]: %s " % (time.time() - start_time))             
## end ##
