import os
from Bio import PDB

def extract_protonation_states(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    protonation_states = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = (chain.id, residue.id[1], residue.resname)
                if residue.resname in ['HIS', 'GLU', 'ASP', 'LYS', 'ARG']:
                    # Check for hydrogen atoms indicating protonation state
                    hydrogen_atoms = [atom for atom in residue if atom.element == 'H']
                    if residue.resname == 'HIS':
                        # Histidine protonation (checking for hydrogen on ND1 and/or NE2)
                        nd1_h = any(atom.name == 'HD1' for atom in hydrogen_atoms)
                        ne2_h = any(atom.name == 'HE2' for atom in hydrogen_atoms)
                        if nd1_h and ne2_h:
                            protonation_states[res_id] = 'HIS (both)'
                        elif nd1_h:
                            protonation_states[res_id] = 'HIS (ND1)'
                        elif ne2_h:
                            protonation_states[res_id] = 'HIS (NE2)'
                        else:
                            protonation_states[res_id] = 'HIS (none)'
                    elif residue.resname == 'GLU':
                        # Glutamate protonation (checking for hydrogen on OE1 and/or OE2)
                        oe1_h = any(atom.name == 'HOE1' for atom in hydrogen_atoms)
                        oe2_h = any(atom.name == 'HOE2' for atom in hydrogen_atoms)
                        if oe1_h or oe2_h:
                            protonation_states[res_id] = 'GLUH'
                        else:
                            protonation_states[res_id] = 'GLU'
                    elif residue.resname == 'ASP':
                        # Aspartate protonation (checking for hydrogen on OD1 and/or OD2)
                        od1_h = any(atom.name == 'HOD1' for atom in hydrogen_atoms)
                        od2_h = any(atom.name == 'HOD2' for atom in hydrogen_atoms)
                        if od1_h or od2_h:
                            protonation_states[res_id] = 'ASPH'
                        else:
                            protonation_states[res_id] = 'ASP'
                    elif residue.resname == 'LYS':
                        # Lysine protonation (checking for hydrogen on NZ)
                        nz_h = any(atom.name == 'HZ3' for atom in hydrogen_atoms)
                        if nz_h:
                            protonation_states[res_id] = 'LYSH'
                        else:
                            protonation_states[res_id] = 'LYS'
                    elif residue.resname == 'ARG':
                        # Arginine protonation (checking for hydrogen on NE, NH1, and NH2)
                        ne_h = any(atom.name == 'HE' for atom in hydrogen_atoms)
                        nh1_h = any(atom.name == 'HH11' for atom in hydrogen_atoms)
                        nh2_h = any(atom.name == 'HH21' for atom in hydrogen_atoms)
                        if ne_h or nh1_h or nh2_h:
                            protonation_states[res_id] = 'ARGH'
                        else:
                            protonation_states[res_id] = 'ARG'
    
    return protonation_states

def compare_protonation_states(pdb_files):
    protonation_states_list = [extract_protonation_states(pdb) for pdb in pdb_files]
    all_residues = set()
    for states in protonation_states_list:
        all_residues.update(states.keys())
    comparison = {res: [states.get(res, '-') for states in protonation_states_list] for res in all_residues}
    return comparison

# Define the directory and file paths
directory = '/home/storm/INBRE-2024/conseneus_protonation_software/'
pdb_files = [os.path.join(directory, 'k1_propka.pqr'), 
             os.path.join(directory, 'k1_pypka.pdb'), 
             os.path.join(directory, 'k1_H++.pdb')]

comparison = compare_protonation_states(pdb_files)

# Formatting and printing the comparison results
header = f"{'Residue':<20} {'k1_propka.pqr':<15} {'k1_pypka.pdb':<15} {'k1_H++.pdb':<15}"
print(header)
print('-' * len(header))
for res, states in comparison.items():
    res_str = f"{res[0]} {res[1]} {res[2]}"
    print(f"{res_str:<20} {states[0]:<15} {states[1]:<15} {states[2]:<15}")
