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

#### Import the starting structure from a file ###############
p = tsase.io.read_con('cluster_2.con')

#### Define the PES ##########################################
lj = tsase.calculators.lj(cutoff=35.0)
p.center(100.0)
p.set_calculator(lj)

###############################################################
#### MOLECULAR DYNAMICS #######################################
###############################################################

#### a few important functions for MD  #####
#### VELOCITY-VERLET ALGORITHM ##############
def step(p,dt,f,fixcm=True):
    m = p.get_momenta()
    m += 0.5 * dt * f

    if fixcm:
                msum = m.sum(axis=0) / float(len(m))
                m = m - msum

    p.set_positions(p.get_positions() + dt * m / p.get_masses()[:,numpy.newaxis])
    p.set_momenta(m)
    f = p.get_forces()
    p.set_momenta(p.get_momenta() + 0.5 * dt * f)
    return f

###################################################################################
##### below is a function which takes a Metropolis Monte Carlo step  ##############
##### The function returns True if the MD step is accepted and False otherwise ####
def MC(p,dr,kT):
    ''' 
    Variables above are as follows
    p = atoms object
    dr = random displacement
    kT = temperature 
    '''
    p_original = p.get_positions()      ## Get original positions
    eng_original = p.get_potential_energy()  ### Get original energy 
    displacement = numpy.random.uniform(-1,1, (len(p.get_positions()),3)) ## Get vector of random displacement
    p_new = p_original + dr*displacement  ### Get new positions
    p.set_positions(p_new)               ## Set new positions
    if p.get_potential_energy() < eng_original:  ### if new energy is less than original accept
        return True
    elif numpy.exp( (eng_original-p.get_potential_energy()) /kT) > numpy.random.uniform(): ## otherwise accept if
        return True
    else:         ## otherwise reject move 
        p.set_positions(p_original)
        return False

################################################################
##### specify MD parameters (Note: In order to convert time into 
##### fs you need to multiple by ase.units.fs

T = 1000.      ### temperature
dt_fs = 0.5    ### time step in fs
dt = dt_fs * ase.units.fs  ## converted time step for md simulations
kT = T * ase.units.kB  

#### this is how to call the thermostat ########################
alpha = 0.8
tcol = 10. * ase.units.fs
therm = tsase.md.nvtandersen(p,dt,kT,alpha,tcol,fixcm=False)

##### Number of MD or MC steps #############################
number_steps = 20000.

##### initiate and run MD trajectory ###########################
f = p.get_forces()
tsase.io.write_con('movie_nvt.con',p,w='w')
timear = []
pe = []


###!!!!!!! PORTION OF CODE TO EDIT BELOW !!!!!!!!!!##########################
#############################################################################
#### create histogram for the number of times a bond lengthed is sampled ####
#############################################################################
tcount=0
fcount=0
mcStepSize=.1
rarray = numpy.zeros(number_steps)
for i in range(int(number_steps)):
    therm.apply_thermostat()
    f = step(p,dt,f)
    #f = MC(p,mcStepSize,kT)
#    if f:
#    	tcount+=1
#    else:
#    	fcount+=1
    if i%100 == 0: ### get snap shot for movie every 100 timesteps 
        tsase.io.write_con('movie_nvt.con',p,w='a')
    r = p.get_positions()
    rarray[i] = numpy.sqrt(numpy.vdot(r[1]-r[0],r[1]-r[0]))

#print tcount, " ", fcount, " ", tcount*1.0/(tcount+fcount)
###!!!!!!!! DO NOT EDIT BELOW THIS LINE !!!!!!!!!!!##########################
#############################################################################
#### Plot the Boltzmann distribution with respect to bond length ############
#### Note: You will not need to change this portion of the code  ############
#############################################################################
numsamp = 100
maxr = 1.5
stepsize = (maxr - 1.0)/numsamp
p.set_positions([[50,50,50],[49.0,50,50]])
r_boltzmann = numpy.zeros(numsamp)
prob_boltzmann = numpy.zeros(numsamp)
norm = 0

for i in range(100):
    x = p.get_positions()
    r_boltzmann[i] = numpy.sqrt(numpy.vdot(x[1]-x[0],x[1]-x[0]))
    prob_boltzmann[i] = numpy.exp(-p.get_potential_energy()/kT)
    if prob_boltzmann[i] > 1e-15:
        norm += stepsize* prob_boltzmann[i]
    x[1][0] -= stepsize
    p.set_positions(x)

prob_boltzmann /= norm

############## Make a plot of the historgram of bondlengths for a diatomic and the probablility distribution of bondlenths from the Boltzmann distribution####### 
figure()
hist(rarray,50,normed=1)
plot(r_boltzmann,prob_boltzmann,c='r')
xlabel('Bond Length (Angstroms)')
ylabel('Probability Density')
savefig('histogram_MD.png')


sys.exit()








