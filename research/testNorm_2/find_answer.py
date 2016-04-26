#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase
from basin import *

##### Setting up the atoms object with ASE and TSASE #########
##############################################################



def basin_hop(p):
    displacement = 0.3
    opt = BasinHopping(  atoms = p, 
                         temperature = 8000 * ase.units.kB,
                         fmax = 0.01,
                         dr = displacement,
                         adjust_cm = True,
                         minenergy = -173.928,
                         #minenergy = -1,
                         distribution = 'uniform',
                         optimizer_logfile = None # uncomment this line if you do not want the information from each optimization to print 

                         )

    opt.record_norm_dists(opt.positions)
    print opt.norm_dists

tcalls=0
for i in range(1): 
    p = tsase.io.read_con('answer_38_less.con')

    lj = tsase.calculators.lj(cutoff=35.0)
    p.center(100.0)
    p.set_calculator(lj)

    basin_hop(p)

    print 'force calls: ', lj.force_calls
    print 'pe: ', p.get_potential_energy()

    tcalls+=lj.force_calls

print "avg: ", tcalls/10

sys.exit()

