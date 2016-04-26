#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase

##### Setting up the atoms object with ASE and TSASE #########
##############################################################

#### Import the starting structure from a file ###############
#p = tsase.io.read_con('diatomic.con')    #lj38-clusters/cluster_0001.con')



def optimize(p):
    p.center(50.0)
    p.set_calculator(lj)

    print '        step   time    potential energy(eV)  force(eV/Angstrom)'
    dmin = tsase.optimize.SDLBFGS(p,maxstep=0.2)
    dmin.run(fmax=0.01,steps=10000) # The number fmax is your convergence criteria and steps is the maxnumber of optimization steps 

    return lj.force_calls - 2

    #tsase.io.write_con('optimized_diatomic_lbfgs.con',p,w='w')

    #print "Pos: ", p.get_positions()
    #print "Dist: ", p.get_all_distances()

curCluster = 0
nArr = []
while curCluster < 100:
    lj = tsase.calculators.lj(cutoff=3.5)
    p = tsase.io.read_con('lj38-clusters/'+str(curCluster)+'.con')
    count = optimize(p)
    nArr.append(count)
    print "PE for ", curCluster, ": ", p.get_potential_energy()
    curCluster+=1

print nArr
maxCount = 0
normTotal = 0
for y in nArr:
    if y >= 10000:
    	maxCount += 1
    else:
    	normTotal += y
print "Hit the max value ", maxCount, " time(s)"
print "Mean of steps for successful tests: ", normTotal / (100-maxCount)


sys.exit()
