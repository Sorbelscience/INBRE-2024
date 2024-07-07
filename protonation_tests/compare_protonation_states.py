# compare_protonation_states.py

import re

def extract_protonation_states(pqr_file):
    protonation_states = {}
    try:
        with open(pqr_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    parts = line.split()
                    atom_name = parts[2]
                    res_name = parts[3]
                    res_id = parts[4]
                    chain_id = parts[5]
                    key = (res_name, res_id, chain_id)
                    if key not in protonation_states:
                        protonation_states[key] = []
                    protonation_states[key].append(atom_name)
        print(f"Extracted protonation states for {pqr_file}")
    except Exception as e:
        print(f"Error reading file {pqr_file}: {e}")
    return protonation_states

def compare_protonation_states(states1, states2):
    differences = {}
    all_keys = set(states1.keys()).union(set(states2.keys()))
    for key in all_keys:
        atoms1 = set(states1.get(key, []))
        atoms2 = set(states2.get(key, []))
        if atoms1 != atoms2:
            differences[key] = (atoms1, atoms2)
    return differences

# Load the protonation states from both PQR files
print("Loading protonation states from PQR files...")
states_4_6 = extract_protonation_states('k1_fixed_ph4.6.pqr')
states_7_0 = extract_protonation_states('k1_fixed_ph7.0.pqr')

# Compare the protonation states
print("Comparing protonation states...")
differences = compare_protonation_states(states_4_6, states_7_0)

# Print the differences
print("Differences in protonation states:")
for key, (atoms_4_6, atoms_7_0) in differences.items():
    res_name, res_id, chain_id = key
    print(f"Residue {res_name} {res_id} in chain {chain_id}:")
    print(f"  pH 4.6: {atoms_4_6}")
    print(f"  pH 7.0: {atoms_7_0}")

if not differences:
    print("No differences found in protonation states between the two pH levels.")
