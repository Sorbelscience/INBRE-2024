from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import argparse

# Argument parser for command line inputs
parser = argparse.ArgumentParser(description='Run OpenMM simulation.')
parser.add_argument('-i', required=True, help='Input CIF file name without extension')
parser.add_argument('-o', required=True, help='Output CIF file name without extension')
args = parser.parse_args()

# Load the system from the CIF file
crd = PDBxFile(args.i + '.cif')

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(crd.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)

# Define the integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Define the simulation
simulation = Simulation(crd.topology, system, integrator)
simulation.context.setPositions(crd.positions)

# Minimize energy
simulation.minimizeEnergy()

# Define the reporters to write outputs
num_snapshots = 750
interval = int(1 * 500000)  # 1 ns intervals, 500000 steps/ns with a 2 fs timestep

# DCD reporter for trajectory
simulation.reporters.append(DCDReporter('trajectory.dcd', interval))

# State data reporter for energy, temperature, and progress
simulation.reporters.append(StateDataReporter(sys.stdout, 10000, step=True, time=True, 
                                              potentialEnergy=True, temperature=True, progress=True, 
                                              remainingTime=True, speed=True, totalSteps=num_steps, separator='\t'))

# Set the number of steps for a 750 ns simulation with a 2 fs timestep
num_steps = int(750 * 1000000 / 2)  # 750 ns * 1,000,000 ns/µs / 2 fs/step

# Run the simulation
simulation.step(num_steps)

# Write the final positions to the CIF file
PDBxFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(args.o + '.cif', 'w'))
