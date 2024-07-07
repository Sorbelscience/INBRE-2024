import parmed as pmd
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# Load the CIF file using ParmEd
structure = pmd.load_file('k1_H++.cif')

# Convert ParmEd structure to OpenMM topology and positions
topology = structure.topology
positions = structure.positions

# Prepare the force field and system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)

# Define the integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)  # Reduced timestep

# Define the simulation
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

# Minimize energy more thoroughly
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=10000)

# Define the reporters to write outputs
num_snapshots = 1000
interval = int((100 * 1000000) / (num_snapshots * 2))  # 100 ns / 1000 snapshots = 0.1 ns; 1 ns = 500000 steps

# Set the number of steps for a 100 ns simulation with a 2 fs timestep
num_steps = int(100 * 1000000 / 2)  # 100 ns * 1,000,000 ns/µs / 2 fs/step

# DCD reporter for trajectory
simulation.reporters.append(DCDReporter('trajectory.dcd', interval))

# State data reporter for energy, temperature, and progress
simulation.reporters.append(StateDataReporter(sys.stdout, interval, step=True, time=True, 
                                              potentialEnergy=True, temperature=True, progress=True, 
                                              remainingTime=True, speed=True, totalSteps=num_steps, separator='\t'))

# Checkpoint reporter
simulation.reporters.append(CheckpointReporter('checkpoint.chk', interval))

# Run the simulation with try-except to catch potential NaN errors
try:
    simulation.step(num_steps)
except Exception as e:
    print(f"Simulation failed with error: {e}")
    print("Check the initial configuration, constraints, and force field parameters.")

# Write the final positions to the CIF file
PDBxFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open('output.cif', 'w'))
