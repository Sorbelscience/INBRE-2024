import re
import os
import shutil
from Bio import PDB
from pdbfixer import PDBFixer
from openmm.app import PDBFile

input_file = "k1_no_H.pdb"
ph_value = 4.6
intermediate_file = "k1_no_H_temp.pdb"
output_pdb = "k1_propka_protonated.pdb"
pka_file = "k1_no_H.pka"

def run_propka(input_pdb, ph):
    command = f"propka3 --pH {ph} {input_pdb}"
    subprocess.run(command, shell=True, check=True)
    print(f"PROPKA3 run completed.")

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
        # Add hydrogens for protonated state
        if residue.resname in ['ASP', 'GLU']:
            residue.resname = 'ASH' if residue.resname == 'ASP' else 'GLH'
        elif residue.resname == 'HIS':
            residue.resname = 'HIP'
        # Other protonation state changes can be added as needed
    elif state == 'deprotonated':
        # Remove hydrogens for deprotonated state
        if residue.resname in ['CYS', 'TYR']:
            residue.resname = 'CYM' if residue.resname == 'CYS' else 'TYM'
        # Other deprotonation state changes can be added as needed

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
    io.save('k1_protonated.pdb')

def add_hydrogens_with_pdbfixer(input_pdb, output_pdb):
    fixer = PDBFixer(filename=input_pdb)
    fixer.addMissingHydrogens(ph=ph_value)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
    print(f"PDBFixer run completed and hydrogens added.")

if __name__ == "__main__":
    input_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", input_file)
    intermediate_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", intermediate_file)
    pka_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", pka_file)
    final_output_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", output_pdb)
    
    # Copy the input file to create an intermediate file
    shutil.copyfile(input_path, intermediate_path)
    
    # Run PROPKA on the intermediate file
    run_propka(intermediate_path, ph_value)
    
    # Parse the PROPKA output to get protonation states
    protonation_states = parse_propka_output(pka_path)
    states_at_pH_4_6 = determine_protonation_states(protonation_states, pH=ph_value)
    
    # Apply protonation states to the PDB file
    apply_protonation_states(intermediate_path, states_at_pH_4_6)
    
    # Add hydrogens using PDBFixer based on modified protonation states
    add_hydrogens_with_pdbfixer('k1_protonated.pdb', final_output_path)
    
    print(f"Final protonated PDB file saved as {final_output_path}")
