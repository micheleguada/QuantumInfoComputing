## Exercise 4 - Quantum Information and Computing 2020/21 ##
# This script change the matrices size in a file and call the FORTRAN program that
# read the file and perfom a matrix-matrix multiplication with three different 
# methods and optimization flags.
#
# AUTHOR: Michele Guadagnini - ID 1230663

import os
import sys
import time
from subprocess import call

#file where the fortran program print fatal messages.
errorfile = "FATAL_ERRORS.txt"  

# flag to activate debug in fortran program
debug = "no"#"TRUE"

sizefile = "MatrixSize.txt"    #file to store the size
flags = ["O0","O2","Ofast"]
folderpath = "Data"
methods = ["1stLoop","2ndLoop","MATMUL"]

Nmin = 50
Nmax = 2001
Incr = 50 #increment

## begin ##
start_time = time.time()
print("Matrix-matrix multiplication test")

# resetting the errorfile
open(errorfile,"w").close()
# ensure existence of the data folder
if not os.path.isdir(folderpath):
    os.makedirs(folderpath)
# resetting the output files:
print("  resetting the results files...", end='')
for jj in range(3): 
    for ii in range(3):
        open(folderpath+"/"+methods[ii]+flags[jj]+".txt", "w").close()
print(" Done.")

for N in range(Nmin,Nmax,Incr):
    # overwriting the size file
    print("\nMatrix size set to "+str(N))
    with open(sizefile, "w") as file1:
        file1.write(str(N))

    # Matrix-matrix multiplications and data storing    
    for jj in range(3):
        call([ "./MatrixTest"+flags[jj]+".x", sizefile, 
              folderpath+"/"+methods[0]+flags[jj]+".txt", 
              folderpath+"/"+methods[1]+flags[jj]+".txt", 
              folderpath+"/"+methods[2]+flags[jj]+".txt", debug ])
        #check for errors
        if (os.stat(errorfile).st_size != 0):
            sys.exit("An error occured during the fortran program execution, causing it to stop. Please check "+errorfile+" for details. Exiting...")
        else:
            print("  multiplication with -"+flags[jj]+" flag successful")
            
print("Time elapsed from start [s]: %s " % (time.time() - start_time)) 
## end ##
#With Nmin=50, Nmax=2001, N=50: time elapsed = 3570.927663564682 seconds
