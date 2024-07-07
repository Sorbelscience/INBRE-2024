import re
import os
from Bio import PDB

input_pdb_file = "k1.pdb"
pka_file = "k1_amber_fixed.pka"  # Relative path to the .pka file
output_pdb_file = "k1_propka_protonated.pdb"
ph_value = 4.6

def parse_propka_output(filename):
    protonation_states = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if re.match(r'^\s*[A-Z]{3}\s+\d+\s+[A-Z]\s+\d+\.\d+', line):
                    parts = line.split()
                    residue_name = parts[0]
                    residue_number = int(parts[1])
                    chain_id = parts[2]
                    pka_value = float(parts[3])
                    protonation_states[(residue_name, chain_id, residue_number)] = pka_value
        print("Protonation states parsed:", protonation_states)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
    return protonation_states

def determine_protonation_states(protonation_states, pH=4.6):
    states = {}
    try:
        for res, pka in protonation_states.items():
            if pka < pH:
                states[res] = 'deprotonated'
            else:
                states[res] = 'protonated'
        print("Determined protonation states:", states)
    except Exception as e:
        print(f"Error determining protonation states: {e}")
    return states

def modify_residue_protonation(residue, state):
    if state == 'protonated':
        if residue.resname in ['ASP', 'GLU']:
            residue.resname = 'ASH' if residue.resname == 'ASP' else 'GLH'
        elif residue.resname == 'HIS':
            # Add both protonated forms of HIS (HD1 and HE2)
            residue.resname = 'HIP'
    elif state == 'deprotonated':
        if residue.resname in ['CYS', 'TYR']:
            residue.resname = 'CYM' if residue.resname == 'CYS' else 'TYM'

def apply_protonation_states(pdb_filename, protonation_states):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('K1', pdb_filename)

    for model in structure:
        for chain in model:
            for residue in chain:
                res_key = (residue.resname, chain.id, residue.id[1])
                if res_key in protonation_states:
                    modify_residue_protonation(residue, protonation_states[res_key])

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('k1_modified.pdb')

if __name__ == "__main__":
    # Paths
    input_path = os.path.join(os.getcwd(), input_pdb_file)
    pka_path = os.path.join(os.getcwd(), pka_file)  # Relative path to the .pka file
    output_path = os.path.join(os.getcwd(), output_pdb_file)

    try:
        # Parse the PROPKA output to get pKa values
        protonation_states = parse_propka_output(pka_path)
        
        # Determine the protonation states at the given pH
        states_at_pH_4_6 = determine_protonation_states(protonation_states, pH=ph_value)
        
        # Apply protonation states to the PDB file
        apply_protonation_states(input_path, states_at_pH_4_6)
        
        # Move the final file to the output path
        os.rename('k1_modified.pdb', output_path)
        
        print(f"Final protonated PDB file saved as {output_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
