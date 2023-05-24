# Python script to run the program regarding Ising Hamiltonian on 
# different values of the parameters: Lambda, N, k, Debug
# It is part of Week 9-10 assignment.
# 
# AUTHOR: Michele Guadagnini - ID 1230663

import os
from glob import glob
import sys
import time
from subprocess import call
import numpy as np
import csv

# parameters:
config     = "Config_Pars.txt"
times_file = "Timings.txt"
errorfile  = "FATAL_ERRORS.txt"
folder     = "Data"

Nset    = np.array([2,3,4,5,6,7,8,9,10,11,12])
Lam_set = np.linspace(0, 3, 31)
kk      = 8
deb     = "F"

# variables for timings
MeanTimes = []

# ensure existence of the data folder
if not os.path.isdir(folder):
    os.makedirs(folder)

## Reset output files
#print("# Python #: Resetting the results files...", end='')
#for ff in glob(folder+"/EigVals*.dat"):
    #os.remove(ff)
#print(" Done.\n")


## begin ##
start_TOTtime = time.time()

# Loop over different Ns
for N in Nset:
    times = []
    
    # Loop over different Lambdas
    for Lam in Lam_set:
        
        # Change the parameters values in config file
        with open(config, "w") as conf:
            conf.write("NN    = "+str(N)  +"\n")
            conf.write("kk    = "+str(kk) +"\n")
            conf.write("Lambda= "+str(Lam)+"\n")
            conf.write("Debug = "+deb     +"\n")
           
        # Run the program
        start_time = time.time()
        call( "./Ising1D.x" )
        times.append(time.time() - start_time)        
        
        #check for errors
        if (os.stat(errorfile).st_size != 0):
            sys.exit("\n# Python #: An error occured during the fortran program execution,"+ 
                    "causing it to stop. Please check "+errorfile+" for details. Exiting...")
        else:
            print("\n# Python #: Computation successful with N="+str(N)+" and Lambda="+str(Lam)+"\n")
    
    MeanTimes.append(np.mean(times))
            
TOTtime = time.time() - start_TOTtime
## end ##

# storing the timings
with open(times_file, "w") as tt:
    tt.write("# N      avg time [s]\n")
    for idx in range(0, len(Nset)):
        tt.write(" %2d      %.6f \n" % (Nset[idx], MeanTimes[idx]) )
    
    tt.write("\n# Total time: " + str(TOTtime) + "\n")

