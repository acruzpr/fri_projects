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
p = tsase.io.read_con('diatomic.con')    #lj38-clusters/cluster_0001.con')

#### Define the PES ##########################################
lj = tsase.calculators.lj(cutoff=3.5)
p.center(50.0)
p.set_calculator(lj)

#### Do a geometry optimization using the FIRE method  #######
print '        step   time    potential energy(eV)  force(eV/Angstrom)'
dmin = tsase.optimize.SDLBFGS(p,maxstep=0.2)
dmin.run(fmax=0.01,steps=10000) # The number fmax is your convergence criteria and steps is the maxnumber of optimization steps 

#### Print the number of force calls required to find the local minimum: This is the total number of force_calls performed in the script
print lj.force_calls

#### Create a new file called 'min.con' with the new optized structure 
tsase.io.write_con('optimized_diatomic_lbfgs.con',p,w='w')

print "Pos: ", p.get_positions()
print "Dist: ", p.get_all_distances()

sys.exit()








