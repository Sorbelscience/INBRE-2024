import re

def parse_propka_output(propka_output_file):
    protonation_states = {}
    with open(propka_output_file, 'r') as f:
        for line in f:
            if re.match(r'^\s*\d+\s+[A-Z]{3}\s+[A-Z]\s+\d+\s+\d+\.\d+\s+\d+\.\d+', line):
                parts = line.split()
                residue_name = parts[1]
                chain_id = parts[2]
                residue_number = int(parts[3])
                pka_value = float(parts[4])
                protonation_states[(residue_name, chain_id, residue_number)] = pka_value
    return protonation_states

def modify_pdb(pdb_file, protonation_states, output_pdb_file):
    with open(pdb_file, 'r') as f, open(output_pdb_file, 'w') as out_f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()
                residue_number = int(line[22:26].strip())
                key = (residue_name, chain_id, residue_number)
                if key in protonation_states:
                    # Modify the line according to the protonation state
                    # This example just prints the modification, you need to handle it based on your requirements
                    print(f'Modifying {residue_name} at {chain_id}{residue_number} with pKa {protonation_states[key]}')
                    # Here you would modify the line based on your specific requirements
                out_f.write(line)
            else:
                out_f.write(line)

# Example usage
protonation_states = parse_propka_output('k1.pka')
modify_pdb('k1.pdb', protonation_states, 'modified_k1.pdb')
