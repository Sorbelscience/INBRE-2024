import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmmtools.testsystems import LennardJonesFluid
from pymbar import MBAR, timeseries

# Simulation settings
pressure = 80*atmospheres
temperature = 120*kelvin
collision_rate = 5/picoseconds
timestep = 2.5*femtoseconds

# Create a Lennard Jones test fluid
sigma = 3.4*angstrom
epsilon = 0.238 * kilocalories_per_mole
fluid = LennardJonesFluid(sigma=sigma, epsilon=epsilon)
[topology, system, positions] = [fluid.topology, fluid.system, fluid.positions]

# Add a barostat
barostat = MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# Retrieve the NonbondedForce
forces = { force.__class__.__name__ : force for force in system.getForces() }
nbforce = forces['NonbondedForce']

# Make two sets of particles
alchemical_particles = set([0])
chemical_particles = set(range(system.getNumParticles())) - alchemical_particles

# Define the energy function for the CustomNonbondedForce
energy_function = 'lambda*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'
energy_function += 'reff_sterics = sigma*(0.5*(1.0-lambda) + (r/sigma)^6)^(1/6);'
energy_function += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
custom_force = CustomNonbondedForce(energy_function)

# Add lambda as a parameter
custom_force.addGlobalParameter('lambda', 1.0)

# Set sigma and epsilon parameters
custom_force.addPerParticleParameter('sigma')
custom_force.addPerParticleParameter('epsilon')
for index in range(system.getNumParticles()):
    [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
    custom_force.addParticle([sigma, epsilon])
    if index in alchemical_particles:
        # Remove the alchemical particle from the existing NonbondedForce
        nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)

# Set the custom force interactions
custom_force.addInteractionGroup(alchemical_particles, chemical_particles)

# Ensure consistent cutoff settings
cutoff_distance = nbforce.getCutoffDistance()
custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
custom_force.setCutoffDistance(cutoff_distance)

system.addForce(custom_force)

# Create an integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Create a simulation
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

# Minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Collect data
nsteps = 2500
niterations = 500
lambdas = np.linspace(1.0, 0.0, 10) # alchemical lambda schedule
nstates = len(lambdas)
u_kln = np.zeros([nstates, nstates, niterations], np.float64)
kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()

for k in range(nstates):
    for iteration in range(niterations):
        print(f'state {k:5d} iteration {iteration:5d} / {niterations:5d}')
        # Set alchemical state
        simulation.context.setParameter('lambda', lambdas[k])
        # Run some dynamics
        simulation.step(nsteps)
        # Compute energies at all alchemical states
        for l in range(nstates):
            simulation.context.setParameter('lambda', lambdas[l])
            u_kln[k, l, iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kT

# Estimate free energy of Lennard-Jones particle insertion
N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
for k in range(nstates):
    [nequil, g, Neff_max] = timeseries.detectEquilibration(u_kln[k, k, :])
    indices = timeseries.subsampleCorrelatedData(u_kln[k, k, :], g=g)
    N_k[k] = len(indices)
    u_kln[k, :, 0:N_k[k]] = u_kln[k, :, indices].T

# Compute free energy differences
mbar = MBAR(u_kln, N_k)
[DeltaF_ij] = mbar.getFreeEnergyDifferences(compute_uncertainty=False)

print("Free energy change to insert a particle =", DeltaF_ij[nstates-1][0])
