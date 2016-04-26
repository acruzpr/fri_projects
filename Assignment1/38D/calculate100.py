#!/usr/bin/env python

import types
import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase

lj = tsase.calculators.lj(cutoff=3.5)

def magnitude(v):
    mag = numpy.sqrt(numpy.vdot(v,v))
    return mag

#def sum(a):
#    if type(a) is numpy.ndarray:
#        ans = 0.0
#        for y in a:
#            ans += sum(y)
#        return ans
#    return abs(a)

def optimize(p):
    maxstep = 10000
    count = 0
    forces = p.get_forces()
    mag = magnitude(forces)
    limitCount = 0
    while mag > 0.01 * p.get_number_of_atoms() and count < maxstep:
        count += 1
        pos = p.get_positions();
        mult = .01
        #mult = 1
        #forceSum = sum(forces)
        forceMax = numpy.amax(forces)
        if (mult*forceMax > .005):
            limitCount += 1
            mult = .005 / forceMax
        #if (mult * forceSum > .2 * p.get_number_of_atoms()):
        #	limitCount += 1
        #	mult = (.2 * p.get_number_of_atoms())/forceSum

        y = 0
        while y < len(forces):
            i = 0
            while i < len(forces[y]):
                pos[y][i] += mult * forces[y][i]
                i+=1
            y+=1
        p.set_positions(pos)
        forces = p.get_forces()
        #p.set_positions(p.get_positions() + .0001 * p.get_forces())
        mag = magnitude(forces)
        #print "Steps: ", count
        #nArr.append(count+1)
        #fArr.append(mag)
        #peArr.append(p.get_potential_energy())
    print "Limited to .01 angstroms ", limitCount, " time(s)"
    print "Took ", count, " steps out of ", maxstep
    return count

nArr = []
fileNumber = 0
while fileNumber < 100:
    p = tsase.io.read_con('lj38-clusters/' + str(fileNumber) + '.con')
    p.center(50.0)
    p.set_calculator(lj)

    print "Starting PE for ", fileNumber, ": ", p.get_potential_energy()
    #print type(p.get_forces())
    count = optimize(p)
    nArr.append(count)

    #print 'Magnitudes:', fArr
    #print "PE's: ", peArr

    print 'potential energy:', p.get_potential_energy()
    #print 'force:', p.get_forces()
    #print "Pos: ", p.get_positions()
    #print "Dist: ", p.get_all_distances()
    fileNumber += 1
    print ""

print nArr
maxCount = 0
normTotal = 0
for y in nArr:
    if y>=10000:
    	maxCount += 1
    else:
    	normTotal += y
print "Hit the max value ", maxCount, " time(s)"
print "Mean of steps for successful tests: ", normTotal / (100-maxCount)

############################################################
## KEEP THIS PART AT THE BOTTOM OF YOUR SCRIPT
#p = tsase.io.write_con('optimized_diatomic.con',p,w='w') 

# The line above will write out a file called optimized_diatomic.con which is the optimized structure 


#figure()
#plot(nArr, fArr)
#xlabel('Nth Step')
#ylabel('Force Magnitude')
#title('Forces over steps')
#yscale('log')
#savefig('LJ38f.png')
#
#figure()
#plot(nArr, peArr)
#xlabel('Nth Step')
#ylabel('Potential Energy')
#title('Potential Energy over steps')
##yscale('log')
#savefig('LJ38pe.png')

sys.exit()
