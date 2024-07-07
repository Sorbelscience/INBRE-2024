import re
import os
from Bio import PDB

input_pka_file = "/home/storm/INBRE-2024/conseneus_protonation_software/#TEMPLATE.pka"
output_pdb_file = "/home/storm/INBRE-2024/conseneus_protonation_software/##TEMPLATE.pdb"

def parse_pka_file(filename):
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

def create_pdb(protonation_states, output_pdb_file):
    with open(output_pdb_file, 'w') as pdb_file:
        for (residue_name, chain_id, residue_number), state in protonation_states.items():
            pdb_file.write(f"ATOM      1  N   {residue_name} {chain_id}{residue_number:4}       0.000   0.000   0.000  1.00  0.00           N\n")
            pdb_file.write(f"ATOM      2  CA  {residue_name} {chain_id}{residue_number:4}       0.000   0.000   0.000  1.00  0.00           C\n")
            pdb_file.write(f"ATOM      3  C   {residue_name} {chain_id}{residue_number:4}       0.000   0.000   0.000  1.00  0.00           C\n")
            pdb_file.write(f"ATOM      4  O   {residue_name} {chain_id}{residue_number:4}       0.000   0.000   0.000  1.00  0.00           O\n")
            if state == 'protonated':
                pdb_file.write(f"ATOM      5  H   {residue_name} {chain_id}{residue_number:4}       0.000   0.000   0.000  1.00  0.00           H\n")
        print(f"PDB file written to {output_pdb_file}")

if __name__ == "__main__":
    #### Parse the .pka file to get pKa values
    protonation_states = parse_pka_file(input_pka_file)
    
    #### Determine the protonation states at the given pH
    states_at_pH_4_6 = determine_protonation_states(protonation_states, pH=4.6)
    
    #### Create a PDB file with the determined protonation states
    create_pdb(states_at_pH_4_6, output_pdb_file)
