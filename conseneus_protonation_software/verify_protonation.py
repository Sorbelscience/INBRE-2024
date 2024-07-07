import os

def parse_pdb(file_path):
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line)
    return atoms

def count_hydrogens(atoms):
    hydrogen_count = 0
    for atom in atoms:
        atom_name = atom[12:16].strip()
        if atom_name.startswith("H"):
            hydrogen_count += 1
    return hydrogen_count

def verify_protonation(original_pdb, protonated_pdb):
    original_atoms = parse_pdb(original_pdb)
    protonated_atoms = parse_pdb(protonated_pdb)
    
    original_hydrogen_count = count_hydrogens(original_atoms)
    protonated_hydrogen_count = count_hydrogens(protonated_atoms)
    
    return original_hydrogen_count, protonated_hydrogen_count

if __name__ == "__main__":
    original_file = "k1_no_H.pdb"
    protonated_file = "k1_propka_protonated.pdb"

    original_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", original_file)
    protonated_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", protonated_file)
    
    original_hydrogen_count, protonated_hydrogen_count = verify_protonation(original_path, protonated_path)
    
    print(f"Hydrogen atoms in original file ({original_file}): {original_hydrogen_count}")
    print(f"Hydrogen atoms in protonated file ({protonated_file}): {protonated_hydrogen_count}")
    
    if protonated_hydrogen_count > original_hydrogen_count:
        print("Verification successful: Hydrogen atoms were added.")
    else:
        print("Verification failed: No hydrogen atoms were added.")
