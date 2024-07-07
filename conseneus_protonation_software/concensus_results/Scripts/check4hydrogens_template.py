pdb_file_path = '/home/storm/INBRE-2024/conseneus_protonation_software/##TEMPLATE'

def check_hydrogens_in_pdb(file_path):
    hydrogen_count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip().startswith('H'):
                hydrogen_count += 1
    return hydrogen_count

#### Check for hydrogen atoms in the PDB file
hydrogen_count = check_hydrogens_in_pdb(pdb_file_path)

if hydrogen_count > 0:
    print(f"Found {hydrogen_count} hydrogen atoms in the file.")
else:
    print("No hydrogen atoms found in the file.")
