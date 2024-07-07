# updated_manual_inspect_pqr.py

def inspect_residue_protonation(pqr_file, residue_name, residue_id):
    protonated_atoms = []
    try:
        with open(pqr_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    parts = line.split()
                    atom_name = parts[2]
                    res_name = parts[3]
                    res_id = parts[4]
                    if res_name == residue_name and res_id == residue_id:
                        protonated_atoms.append(atom_name)
        return protonated_atoms
    except Exception as e:
        print(f"Error reading file {pqr_file}: {e}")
        return []

# List of residues to inspect
residues_to_check = [
    ("ASP", "57"),
    ("GLU", "54"),
    ("HIS", "126"),
    ("LYS", "56"),
    ("ARG", "62")
]

# Check each residue in both PQR files
for residue_name, residue_id in residues_to_check:
    print(f"Checking {residue_name} {residue_id} at pH 4.6...")
    ph46_protonated_atoms = inspect_residue_protonation('k1_fixed_ph4.6.pqr', residue_name, residue_id)
    print(f"Protonated atoms at pH 4.6: {ph46_protonated_atoms}")
    print(f"Checking {residue_name} {residue_id} at pH 7.0...")
    ph70_protonated_atoms = inspect_residue_protonation('k1_fixed_ph7.0.pqr', residue_name, residue_id)
    print(f"Protonated atoms at pH 7.0: {ph70_protonated_atoms}")
    print("\n")
