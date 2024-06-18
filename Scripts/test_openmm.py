import openmm
from openmm import app
from openmm import unit
from sys import stdout

# Print the OpenMM version
print("OpenMM version:", openmm.version.version)

# Create a simple water box simulation (you need an input PDB file for this example)
pdb = app.PDBFile('input.pdb')
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

# Use NoCutoff method instead of PME
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
