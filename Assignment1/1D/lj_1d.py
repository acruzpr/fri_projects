#!/usr/bin/env python

### importing important libraries 
import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *

###### Lennard-Jones Potential ################
###### The function below returns the potential energy for a given bond length (r), sigma and epsilon
def LJ(r,sigma,epsilon):
    Vlj = 4* epsilon * ( numpy.power(sigma/r,12) - numpy.power(sigma/r,6) )
    return Vlj

###### Defining parameters for your L-J potential
sigma = 1.1
epsilon = 1.0

##### create an array of various bond lengths (rarray) 
##### and its corresponding energy (energyarray)
rmin = 0.9
rmax = 3.0
numbersamples = 100.
rarray = numpy.linspace(rmin,rmax,numbersamples)
energyarray = numpy.zeros(numbersamples)
for i in range(len(rarray)):
    energyarray[i] = LJ(rarray[i],sigma,epsilon)

# The code below plots how the potential energy changes with bondlength using matplotlib 
# Try commenting out the plot and scatter functions below to see how they are different! 

figure()
plot(rarray,energyarray)     #  This will plot a line connecting all points
scatter(rarray,energyarray)  #  This will plot the points sampled
xlabel('Bond Length')
ylabel('Potential Energy')
title('Lennard-Jones potential for diatomic: Sigma=1.1, Epsilon=1')
savefig('LJ_1.1.png')



