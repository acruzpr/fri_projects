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
p = tsase.io.read_con('cluster_38.con')

#### Define the PES ##########################################
lj = tsase.calculators.lj(cutoff=35.0)
p.center(50.0)
p.set_calculator(lj)

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

T = 2000.                     # Set the initial temperature
dt_fs = 1.5                  # Set the initial time step for the MD simulations
dt = dt_fs * ase.units.fs    # Convert time step to correct units
kT = T * ase.units.kB        # Convert temperature to units of energy 

#### this is how to call the thermostat ########################
alpha = 0.8                  # Set Alpha for MD simulation
tcol = 100. * ase.units.fs   # Set tcol for MD simulations
therm = tsase.md.nvtandersen(p,dt,kT,alpha,tcol,fixcm=False)  # Call the thermostat from MD simulation

##### time of MD trajecotory in fs #############################
time = 100000./dt_fs          # Set total time of the MD trajectory (Note: This trajectory is for 50,000 femtoseconds and will work for any timestep)

##### initiate and run MD trajectory ###########################
f = p.get_forces()
tsase.io.write_con('movie_nvt.con',p,w='w')  # initiates a movie for the MD trajectory 

initialEnergy = p.get_total_energy()
energyDiffArr = [0]
bondLenArr = [bondLength(p.get_positions())]
forceArr = [mag(p.get_forces())]
stepNumberArr = [0]
tempArr = [p.get_temperature()]
peArr = [p.get_potential_energy()]

bestP = 100000000

print dt

for i in range(int(time)):
    therm.apply_thermostat()                 # apply the thermostat at each MD step
    f = step(p,dt,f)                         # take an MD step 

    energyDiffArr.append(abs(initialEnergy - (p.get_total_energy())))
    currTime = dt_fs*i
    stepNumberArr.append(currTime)
    bondLenArr.append(bondLength(p.get_positions()))
    forceArr.append(mag(p.get_forces()))
    currP = p.get_potential_energy()
    peArr.append(currP)
    if (currP<bestP):
    	bestP = currP
        print 'Found better P: ', bestP
        tsase.io.write_con('bestP_movie.con',p,w='a')  # initiates a movie for the MD trajectory 
    tempArr.append(p.get_temperature())

    T = abs(1200*(1+math.sin(currTime*50/time))*(1-(currTime/(time*dt_fs))))
    kT = T * ase.units.kB        # Convert temperature to units of energy 
    therm = tsase.md.nvtandersen(p,dt,kT,alpha,tcol,fixcm=False)  # Call the thermostat from MD simulation

    if i%100 == 0:                           ### get snap shot for movie every 100 timesteps 
        tsase.io.write_con('movie_nvt.con',p,w='a')   # Append to movie of the MD trajectory

dmin = tsase.optimize.SDLBFGS(p,maxstep=0.2)
dmin.run(fmax=0.01,steps=10000)

print "Optimized P: ", p.get_potential_energy();
tsase.io.write_con('bestP.con',p,w='w')

print 'done'

figure()
plot(stepNumberArr, energyDiffArr)
xlabel('Time')
ylabel('Energy Deviation')
title('Energy Deviation over time, Timestep=' + str(dt_fs))
savefig('deviation.png')


figure()
plot(stepNumberArr, bondLenArr)
xlabel('Time')
ylabel('Bond Length')
title('Bond Length over time')
savefig('bond_len_deviation.png')

figure()
plot(stepNumberArr, forceArr)
xlabel('Time')
ylabel('Force')
title('Force over time')
savefig('force_deviation.png')

figure()
plot(stepNumberArr, tempArr)
xlabel('Time')
ylabel('Temperature')
title('Temperature over time')
savefig('temp_deviation.png')

figure()
plot(stepNumberArr, peArr)
xlabel('Time')
ylabel('Potential Energy')
title('Potential Energy over time, ' + str(T) + 'K')
savefig('p_deviation.png')


sys.exit()








