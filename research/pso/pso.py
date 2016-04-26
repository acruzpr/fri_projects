#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
from math import *
import tsase
import ase
from basin import BasinHopping
import scipy
import time
from bcm import BCM

#Constants
POPULATION = 4
NUM_BESTS = 3
MIN_ENERGY = 173.928427
BEST_BCM_CUTOFF = -1


def mag(cls,vec):
    return sqrt(dot(vec,vec))

def get_basin_hop(p,energy_goal):
    displacement = 0.3
    opt = BasinHopping(  atoms = p, 
                         temperature = 8000 * ase.units.kB,
                         fmax = 0.01,
                         dr = displacement,
                         minenergy = energy_goal,
                         optimizer_logfile = None # uncomment this line if you do not want the information from each optimization to print 

                         )
    return opt

lj = tsase.calculators.lj(cutoff=35.0)
def random_structure():
    r = tsase.io.read_con(random_structure_arr[int(random()*100)])
    r.center(100.0)
    r.set_calculator(lj)
    ret = get_basin_hop(r,MIN_ENERGY)
    return ret

class placeholder_structure:
    def get_bcm(self):
        return [0] * 6
    def get_energy(self):
        return 0

#inserts into the list of bests if better
def insert_best(bests, structure):
    nenergy = structure.get_energy()
    print 'nenergy: ', nenergy
    if (bests[0].get_energy()>nenergy):
        i=0
        a=0
        within_cutoff=None
        while (a < len(bests)):
            if bests[i].get_energy()>nenergy:
                i = a
            temp_mag = structure.get_bcm()-bests[a].get_bcm()
            if (within_cutoff is None and mag(temp_mag, temp_mag) < BEST_BCM_CUTOFF):
                print 'within_cutoff: ', temp_mag
                within_cutoff = a
            a+=1
            if(a<len(bests)):
                print 'energy: ', bests[a].get_energy()
        print i, within_cutoff
        if (within_cutoff is not None):
            if (bests[within_cutoff].get_energy()>nenergy):
                scopy = structure.copy()
                bests[within_cutoff] = scopy
                while (within_cutoff < i and within_cutoff < len(bests)-1):
                    within_cutoff += 1
                    bests[within_cutoff-1] = bests[within_cutoff]
                    bests[within_cutoff] = scopy
        else:
            temp2 = structure.copy()
            while i>0:
                i-=1
                temp = bests[i]
                bests[i] = temp2
                temp2 = temp


print 'start'
random_structure_arr = [0]*100
for i in range(100):
    random_structure_arr[i] = '../100_lj38/'+str(i)+'.con'

structures = [None] * POPULATION
for i in range(POPULATION):
    structures[i] = random_structure()
    print i, ': ',  structures[i].get_bcm()


bests = [placeholder_structure()] * NUM_BESTS
for i in range(POPULATION):
    insert_best(bests, structures[i])

for i in range(NUM_BESTS):
    print i, ': ', bests[i].get_energy(), " bcm: ", bests[i].get_bcm()



#sys.exit()


