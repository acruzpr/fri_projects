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

#### The line below imports the starting structure from a file called diatomic.con ###############
p = tsase.io.read_con('diatomic.con')

## The lines below set the potential energy surface to be a lennard-jones cluster  ###############
lj = tsase.calculators.lj(cutoff=3.5)
p.center(50.0)
p.set_calculator(lj)

### In this script the variable p is your 'atoms object' or the variable which holds the information about the atomic system #############

def magnitude(v):
    print "mag: ", v
    mag = numpy.sqrt(numpy.vdot(v,v))
    return mag

pos = p.get_positions()
pos[0][0] = 50.95
p.set_positions(pos)

deltax = 0.0 * p.get_forces()
deltax[0][0] = 0.001
print deltax

maxstep = 50000
count = 0
while magnitude(p.get_forces()) > 0.1 and count < maxstep:
	p.set_positions(p.get_positions() + deltax)
	pos = p.get_positions()
	print magnitude(pos[1]-pos[0])
	count += 1

### To get the positions of your molecular system you can type the following #### 
print 'Note: the positions below are the x,y, and z coordinates for each atom'
print 'positions:',p.get_positions()


### To get the potential energy your can type the following  ##### 
print 'potential energy:',p.get_potential_energy()

### To get the forces type the following ###
print 'force:',p.get_forces()





sys.exit()








