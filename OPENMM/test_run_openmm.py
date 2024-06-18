from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Load the PDB file
pdb = PDBFile('1AKI.pdb')

# Load the force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create the system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Set up the integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create the simulation object
simulation = Simulation(pdb.topology, system, integrator)

# Set the initial positions
simulation.context.setPositions(pdb.positions)

# Minimize the energy
simulation.minimizeEnergy()

# Set up reporters to output data
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('output.pdb', 1000))

# Run the simulation for 10,000 steps
simulation.step(10000)

