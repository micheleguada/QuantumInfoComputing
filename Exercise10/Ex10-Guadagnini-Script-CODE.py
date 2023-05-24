# Python script to run the program regarding Renormalization Group
# applied to 1D Ising Model over different values of the parameters:
#   Lambda, N, TruncN ...
# It is part of Week 11-12 assignment (Ex10).
# 
# AUTHOR: Michele Guadagnini - ID 1230663

import os
from glob import glob
import sys
import time
import subprocess
import numpy as np

# python parameters
verbosity = False
folder    = "Data"

# fortran parameters:
config    = "Config_Pars.txt"
errorfile = "FATAL_ERRORS.txt"

Nset    = np.array([4,4,4,6,6,6,8,8,8])
TruncN  = np.array([2,3,4,2,3,4,2,3,4])
Lam_set = np.linspace(0, 3, 61)
ExitTh  = "5d-16"
MaxItr  = 100
kk      = 1
deb     = "F"

# variables for timings
times_file = "Timings.txt"
MeanTimes  = []
StdDev     = []

# ensure existence of the data folder
if not os.path.isdir(folder):
    os.makedirs(folder)

# Reset output files
print("# Python #: Resetting the results files...", end='')
for ff in glob(folder+"/EigVals*.dat"):
    os.remove(ff)
print(" Done.\n")


## begin ##
start_TOTtime = time.time()

# Loop over different Ns
for N,T in zip(Nset,TruncN):
    times = []
    
    # Loop over different Lambdas
    for Lam in Lam_set:
        
        # Change the parameters values in config file
        with open(config, "w") as conf:
            conf.write("NN    = "+str(N)     +"\n")
            conf.write("TruncN= "+str(T)     +"\n")
            conf.write("Lambda= "+str(Lam)   +"\n")
            conf.write("ExitTh= "+ExitTh     +"\n")
            conf.write("MaxItr= "+str(MaxItr)+"\n")
            conf.write("kk    = "+str(kk)    +"\n")
            conf.write("Debug = "+deb        +"\n")
           
        # Run the program
        start_time = time.time()
        if verbosity:
            subprocess.call( "./RGIsing1D.x" )
        else:
            subprocess.run( ["./RGIsing1D.x"], stdout=subprocess.DEVNULL )
        times.append(time.time() - start_time)        
        
        #check for errors
        if (os.stat(errorfile).st_size != 0):
            sys.exit("\n# Python #: An error occured during the fortran program execution,"+ 
                    "causing it to stop. Please check "+errorfile+" for details. Exiting...")
        else:
            print("# Python #: Computation successful with N="+str(N)+", TruncN="+str(T)+" and Lambda="+str(Lam)+"\n")
    
    MeanTimes.append(np.mean(times))
    StdDev.append(np.std(times))
            
TOTtime = time.time() - start_TOTtime
## end ##

# storing the timings
with open(times_file, "w") as tt:
    tt.write("# N   TruncN   avg time [s]  std [s]\n")
    for idx in range(0, len(Nset)):
        tt.write(" %2d      %2d      %.6f      %.6f \n" % (Nset[idx], TruncN[idx], MeanTimes[idx], StdDev[idx]) )
    
    tt.write("\n# Total time: " + str(TOTtime) + "\n")

