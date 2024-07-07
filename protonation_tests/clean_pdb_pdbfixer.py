from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='k1_with_h.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)  # You can change the pH value here if needed
PDBFile.writeFile(fixer.topology, fixer.positions, open('k1_fixed.pdb', 'w'))
