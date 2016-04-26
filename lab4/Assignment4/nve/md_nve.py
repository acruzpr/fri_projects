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

#### Do a geometry optimization using the FIRE method  #######
dmin = ase.optimize.FIRE(p)
dmin.run(fmax=0.01)

###############################################################
###############################################################
#### MOLECULAR DYNAMICS #######################################
###############################################################

def mag(v):
    return numpy.sqrt(numpy.vdot(v,v))

def bondLength(positions):
    return mag(positions[0]-positions[1])

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

################################################################
##### specify MD parameters (Note: In order to convert time into 
##### fs you need to multiple by ase.units.fs

T = 300.                    # Set the initial temperature
dt_fs = 1                 # Set the initial time step for the MD simulations
dt = dt_fs * ase.units.fs   # Convert time step to correct units
kT = T * ase.units.kB       # Convert temperature to units of energy 

################################################################
#### get initial kinetic energy from the boltzmann distribution#
masses = p.get_masses()
vel = [[-numpy.sqrt(kT/masses[0]),0,0],[numpy.sqrt(kT/masses[0]),0,0]] # Set initial velocity; gives diatomic vibrational kinetic energy
p.set_velocities(vel)

##### time of MD trajecotory in fs #############################
time = 20./dt_fs            # Set total time of the MD trajectory (Note: This trajectory is for 20 femtoseconds and will work for any timestep)

##### initiate and run MD trajectory ###########################
f = p.get_forces()
tsase.io.write_con('movie_nve.con',p,w='w') # initiates a movie for the MD trajectory


initialEnergy = p.get_total_energy()
energyDiffArr = [0]
bondLenArr = [bondLength(p.get_positions())]
forceArr = [mag(p.get_forces())]
stepNumberArr = [0]
print 'pos: ', p.get_positions()

for i in range(int(time)):
    f = step(p,dt,f)        # take an MD step 
    #print i
    energyDiffArr.append(abs(initialEnergy - (p.get_total_energy())))
    stepNumberArr.append(dt*(i+1))
    bondLenArr.append(bondLength(p.get_positions()))
    forceArr.append(mag(p.get_forces()))
    if i%1 == 0:
        tsase.io.write_con('movie_nve.con',p,w='a')  # Append to movie of the MD trajectory

print energyDiffArr
print stepNumberArr

figure()
plot(stepNumberArr, energyDiffArr)
xlabel('Time')
ylabel('Energy Deviation')
title('Energy Deviation over time, NVE, timestep=' + str(dt_fs) )
savefig('deviation.png')


figure()
plot(stepNumberArr, bondLenArr)
xlabel('Time')
ylabel('Bond Length')
title('Bond Length over time, NVE')
savefig('bond_len_deviation.png')

figure()
plot(stepNumberArr, forceArr)
xlabel('Time')
ylabel('Force')
title('Force over time, NVE')
savefig('force_deviation.png')

figure()
plot(stepNumberArr, forceArr)
xlabel('Time')
ylabel('Temperature')
title('Temperature over time, NVE')
savefig('temp_deviation.png')

sys.exit()
