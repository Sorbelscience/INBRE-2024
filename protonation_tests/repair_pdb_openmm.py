from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

pdb = PDBFile('cleaned_k1.pdb')
forcefield = ForceField('amber99sb.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=4.6)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print('Minimizing...')
simulation.minimizeEnergy()

positions = simulation.context.getState(getPositions=True).getPositions()
with open('repaired_k1_openmm.pdb', 'w') as output:
    PDBFile.writeFile(modeller.topology, positions, output)

print('Minimization done. Repaired PDB file saved as repaired_k1_openmm.pdb')
