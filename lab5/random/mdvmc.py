#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase

rarray = []

#seed = 0;
#tval=0
#tdone=0
#for i in range(1000):
#    print seed
#    tval+=seed
#    tdone+=1
#    rarray.append(seed)
#    seed = (seed*seed+1)%(65536-seed)

#seed = 123456789
#tval = 0
#tdone = 0
#a = 1103515245
#c = 12345
#m = math.pow(2,32)
#for i in range(1000):
#	print seed
#	tval += seed
#	tdone += 1
#	rarray.append(seed)
#	seed = (seed * a + c) % m


tval = 0
tdone = 0
for i in range(100000):
    num = np.random.rand()
    print num
    ans = math.sqrt(num)/4
    if(0==(int)(num*33554432)&1):
    	ans+=.5
    if(0==(int)(num*16777216)&1):
    	ans = 1-ans
    rarray.append(ans)
    tval += ans
    tdone += 1
    print 'ans: ', ans

print 'avg: ', tval/tdone

figure()
hist(rarray,50,normed=1)
xlabel('Bond Length (Angstroms)')
ylabel('Probability Density')
savefig('histogram_MD.png')


sys.exit()








