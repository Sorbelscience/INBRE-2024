from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer('k1.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens()
PDBFile.writeFile(fixer.topology, fixer.positions, open('k1_fixed.pdb', 'w'))
