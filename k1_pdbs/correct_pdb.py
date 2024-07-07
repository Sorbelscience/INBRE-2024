from Bio.PDB import PDBParser, PDBIO

# Load the PDB files
parser = PDBParser(QUIET=True)
structure_k1 = parser.get_structure('k1', '/home/storm/INBRE-2024/k1_pdbs/k1.pdb')
structure_matk1_1 = parser.get_structure('matk1_1', '/home/storm/INBRE-2024/k1_pdbs/matk1_1.pdb')
structure_matk1_2 = parser.get_structure('matk1_2', '/home/storm/INBRE-2024/k1_pdbs/matk1_2.pdb')

# Function to correct atom names based on reference structures
def correct_atom_names(target_structure, reference_structures):
    for ref_structure in reference_structures:
        for ref_model, target_model in zip(ref_structure, target_structure):
            for ref_chain, target_chain in zip(ref_model, target_model):
                for ref_residue, target_residue in zip(ref_chain, target_chain):
                    ref_atoms = {atom.get_id(): atom for atom in ref_residue}
                    for target_atom in target_residue:
                        if target_atom.get_id() in ref_atoms:
                            ref_atom = ref_atoms[target_atom.get_id()]
                            target_atom.name = ref_atom.name
                            target_atom.element = ref_atom.element

# Apply corrections based on reference structures
correct_atom_names(structure_k1, [structure_matk1_1, structure_matk1_2])

# Function to add TER records
def add_ter_records(structure):
    for model in structure:
        for chain in model:
            last_residue_id = None
            for residue in chain:
                if last_residue_id is not None and residue.get_id()[1] != last_residue_id + 1:
                    chain.child_list.append((' ', 'TER', ' '))
                last_residue_id = residue.get_id()[1]
    return structure

# Add TER records
structure_k1 = add_ter_records(structure_k1)

# Save the corrected structure
io = PDBIO()
io.set_structure(structure_k1)
corrected_pdb_path = '/home/storm/INBRE-2024/k1_pdbs/k1_corrected_renamed.pdb'
io.save(corrected_pdb_path)

print(f'Corrected PDB file saved at {corrected_pdb_path}')
