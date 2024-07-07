from Bio import PDB

def parse_protonation_states(filename):
    protonation_states = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split(', Protonation State: ')
            res_info = eval(parts[0].replace('Residue: ', ''))
            state = parts[1]
            protonation_states[res_info] = state
    return protonation_states

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

# Example usage
protonation_states = parse_protonation_states('protonation_states_pH_4.6.txt')
apply_protonation_states('k1.pdb', protonation_states)

print("Protonated PDB file saved as 'k1_protonated.pdb'")
