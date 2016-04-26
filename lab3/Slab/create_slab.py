#!/usr/bin/env python

import numpy
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import sys
import ase

###################################################################
#### PARAMETERS TO CHANGE ########################################
###################################################################

xyz_repeat = [1.,1.,4.] # number of unit cell repeats in the x, y, and z direction respectively
zaddition = 15.   # additional space vacuum in angstroms
AO = True  # Indicates where you cleave the perovskite surface Either at SrO (True) or FeO2 (False)

###################################################################
## Note Bottom two layers of slab are frozen#######################
###################################################################

f = open('CONTCAR','r')
lines = f.readlines()
f.close()

x = lines[6].split()
numatoms = 0 
numorig = 0
for i in range(len(x)):
    numatoms += float(x[i])*8
    numorig += float(x[i])
origcell_length = float(lines[2].split()[0])



nf = open('POSCAR','w')

for i in range(len(lines)):
    if i < 2:
        nf.write(lines[i])
    elif i < 5:
        
        split = lines[i].split()
        for j in range(len(split)):
         if (j == 2) and (float(split[j]) > 0.1):
            split[j] = xyz_repeat[j] * float(split[j]) + zaddition 
            split[j] = str('%9.16f' % split[j])
         else:
            split[j] = xyz_repeat[j] * float(split[j])
            split[j] = str('%9.16f' % split[j])

        nf.write("    ")
        nf.write(" ".join(split))
        nf.write("\n")
        
    elif i < 6:
        nf.write(lines[i])
    elif i < 7: 
        split = lines[i].split()
        for j in range(len(split)):
            split[j] = numpy.prod(xyz_repeat) * float(split[j])
            split[j] = str(int(split[j]))
        nf.write("   ")
        nf.write(" ".join(split))
        nf.write("\n")
    elif i < 8:
        nf.write("Selective dynamics \n")
        nf.write(lines[i])
    elif i < (8+ numorig):
        zshift = float(xyz_repeat[2]*origcell_length)/(xyz_repeat[2]*origcell_length + zaddition)
        split1 = lines[i].split()
        split = lines[i].split()
        for k in range(int(xyz_repeat[0])):
            for w in range(int(xyz_repeat[1])):
                for q in range(int(xyz_repeat[2])):
                    split1[0] = str('%9.16f' % (1./xyz_repeat[0] * float(split[0]) + (float(k))/xyz_repeat[0]))
                    split1[1] = str('%9.16f' % (1./xyz_repeat[1] * float(split[1]) + (float(w))/xyz_repeat[1]))
                    split1[2] = str('%9.16f' % (1./xyz_repeat[2] * float(split[2]) + (float(q))/xyz_repeat[2]))
                    if AO == False:
                        split1[2] = str('%9.16f' % (zshift*float(split1[2])+(1-zshift)))
                        if float(split1[2]) < (1. - zshift/2.):
                            a = '   T T T'
                        else:
                            a = '   F F F'
                    else:
                        split1[2] = str('%9.16f' % (zshift*float(split1[2])))
                        if float(split1[2]) > (zshift/2. - 0.01):
                            a = '   T T T'
                        else:
                            a = '   F F F'
                    nf.write("   ")
                    nf.write(" ".join(split1))
                    nf.write(a +"\n")

    elif i == (8+ numorig + 1):
        nf.write("\n")
        for j in range(int(numatoms)):
            nf.write(lines[i])

nf.close()

sys.exit()

