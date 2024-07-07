import re

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

def save_protonation_states(states, output_file):
    try:
        with open(output_file, 'w') as f:
            for res, state in states.items():
                f.write(f"Residue: {res}, Protonation State: {state}\n")
        print(f"Protonation states saved to {output_file}")
    except Exception as e:
        print(f"Error writing to file {output_file}: {e}")

# Example usage
protonation_states = parse_propka_output('k1.pka')
if protonation_states:
    states_at_pH_4_6 = determine_protonation_states(protonation_states, pH=4.6)
    save_protonation_states(states_at_pH_4_6, 'protonation_states_pH_4.6.txt')

print(f"Protonation states have been saved to 'protonation_states_pH_4.6.txt'.")
