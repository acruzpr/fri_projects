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

def F(r,sigma,epsilon):
    f = -24*epsilon/r * (numpy.power(sigma/r,6) - 2*numpy.power(sigma/r,12))
    return f


###### Defining parameters for your L-J potential
sigma = 1.0
epsilon = 1.0

##### create an array of various bond lengths (rarray) 
##### and its corresponding energy (energyarray)
#rmin = 0.9
#rmax = 3.0
#numbersamples = 100.
#rarray = numpy.linspace(rmin,rmax,numbersamples)
#energyarray = numpy.zeros(numbersamples)
#for i in range(len(rarray)):
#    energyarray[i] = LJ(rarray[i],sigma,epsilon)

rAt = .95
a = .0001
goal = .001
n = 1
nArr = [n]
rArr = [rAt]
magArr = [F(rAt, sigma, epsilon)]
peArr = [LJ(rAt, sigma, epsilon)]

f = F(rAt, sigma, epsilon)
lj = 0
while f > goal:
    # print "R: %f, F: %f, LJ: %f" % (rAt, f, lj)
    rAt = rAt + a * f
    n = n + 1 
    nArr = numpy.append(nArr, n)
    rArr = numpy.append(rArr, rAt)
    f = F(rAt,sigma,epsilon)
    magArr = numpy.append(magArr, f)
    lj = LJ(rAt,sigma,epsilon)
    peArr = numpy.append(peArr, lj)

print rAt
print nArr

# The code below plots how the potential energy changes with bondlength using matplotlib 
# Try commenting out the plot and scatter functions below to see how they are different! 

figure()
plot(nArr, rArr)     #  This will plot a line connecting all points
xlabel('Nth Step')
ylabel('Radius')
title('Radius used in iterative approximation, a=.00001')
savefig('LJ7R.png')

figure()
plot(nArr, magArr)     #  This will plot a line connecting all points
xlabel('Nth Step')
ylabel('Force')
title('Forces in iterative approximation, a=.00001')
savefig('LJ7F.png')

figure()
plot(nArr, peArr)     #  This will plot a line connecting all points
xlabel('Nth Step')
ylabel('Potential Energy')
title('Potential Energy in iterative approximation, a=.00001')
savefig('LJ7E.png')
