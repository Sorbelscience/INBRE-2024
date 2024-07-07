from openmm.app import PDBFile
from openmm.unit import *

def get_pdb_dimensions(pdb_file):
    pdb = PDBFile(pdb_file)
    positions = pdb.positions
    min_x, min_y, min_z = float('inf'), float('inf'), float('inf')
    max_x, max_y, max_z = float('-inf'), float('-inf'), float('-inf')

    for pos in positions:
        if pos.x < min_x:
            min_x = pos.x
        if pos.y < min_y:
            min_y = pos.y
        if pos.z < min_z:
            min_z = pos.z
        if pos.x > max_x:
            max_x = pos.x
        if pos.y > max_y:
            max_y = pos.y
        if pos.z > max_z:
            max_z = pos.z

    dimensions = (max_x - min_x, max_y - min_y, max_z - min_z)
    return dimensions

# Use the function to get dimensions of your PDB file
pdb_file = 'k1_H++.pdb'
dimensions = get_pdb_dimensions(pdb_file)
print(f"Dimensions of the PDB file: {dimensions}")
