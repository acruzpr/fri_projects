#!/usr/bin/env python

import numpy
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import sys

##################################################################
## Replace with your list of energy and lattice constant here ####
##################################################################


latticec = [3.75,3.8,3.85,3.9,3.95,4,4.05,4.1]
energy = [31.039,31.279,31.412,31.468,31.452,31.377,31.238,31.049]

##################################################################
##################################################################
##################################################################

info = numpy.polyfit(latticec,energy,2)  # find quadratic fit to your data

### plot your data and your quadratic line of best fit ##########
lc = numpy.linspace(latticec[0],latticec[len(latticec)-1],100)
quadfit = numpy.zeros(len(lc))
for i in range(len(lc)):
    quadfit[i] = numpy.polyval(info,lc[i])
figure()
plot(lc,quadfit,c='r')
scatter(latticec,energy,c='k')
xlabel('lattice constant (Angstroms)')
ylabel('potential energy (eV)')
savefig('latticefit.png')
#################################################################

### Find the minimum by taking the derivative (polyder) and finding its zeros (roots) ####
der = numpy.roots(numpy.polyder(info))
print 'your optimal lattice constant:',der[0]
sys.exit()

