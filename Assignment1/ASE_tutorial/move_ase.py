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

print 'old positions:', p.get_positions()
new_positions = p.get_positions() + 0.1* p.get_forces()

###########################################################
# To move the positions do the following: #################

p.set_positions(new_positions)
print 'new positions:', p.get_positions()





sys.exit()








