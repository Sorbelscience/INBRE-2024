from Bio.PDB import PDBParser

# Load the PDB files
parser = PDBParser(QUIET=True)
structure_input_fix = parser.get_structure('input_fix', '/home/storm/INBRE-2024/k1_pdbs/input.fix.pdb')
structure_k1 = parser.get_structure('k1', '/home/storm/INBRE-2024/k1_pdbs/k1.pdb')

# Function to extract atom information
def extract_atom_info(structure):
    atom_info = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord_tuple = tuple(atom.coord)
                    atom_info.add((chain.id, residue.get_resname(), residue.id, atom.get_id(), coord_tuple))
    return atom_info

# Extract atom information from both structures
atoms_input_fix = extract_atom_info(structure_input_fix)
atoms_k1 = extract_atom_info(structure_k1)

# Find unique atoms
unique_to_input_fix = atoms_input_fix - atoms_k1
unique_to_k1 = atoms_k1 - atoms_input_fix

# Print unique atoms
print("Unique to input.fix.pdb:")
for atom in unique_to_input_fix:
    print(atom)

print("\nUnique to k1.pdb:")
for atom in unique_to_k1:
    print(atom)
