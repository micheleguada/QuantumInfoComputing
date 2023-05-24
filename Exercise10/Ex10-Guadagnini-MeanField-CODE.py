# Python script to compute the mean field solution of
# the ground state in function of Lambda for the 1D
# Ising model with transverse field.
# It is part of Week 11-12 assignment (Ex10).
# 
# AUTHOR: Michele Guadagnini - ID 1230663

import os
import numpy as np

LambdaSet   = np.linspace(0, 3, 61)
MeanFieldGS = "Data/MeanFieldGS.dat"
folder      = "Data"

def MeanField(Lam):
    if ((Lam > -2) and (Lam < 2)):
        ee = -1 - Lam**2/4
    else:
        ee = - abs(Lam)
    return ee
   
# ensure existence of the data folder
if not os.path.isdir(folder):
    os.makedirs(folder)  

## begin ##    
# saving the mean field
with open(MeanFieldGS, "w") as mf:
    mf.write("# Lam   GS_eigval \n")
    for Lam in LambdaSet:
        mf.write("%.8f      %.8f \n" % ( Lam, MeanField(Lam) ))

## end ##
