import os

def parse_pdb(file_path):
    #### Parses a PDB file to extract lines that start with "ATOM" or "HETATM".
    #
    # Parameters:
    # file_path (str): The path to the PDB file to be parsed.
    #
    # Returns:
    # list: A list of lines from the PDB file that represent atoms.
    ####
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line)
    return atoms

def count_hydrogens(atoms):
    #### Counts the number of hydrogen atoms in a list of atom lines from a PDB file.
    #
    # Parameters:
    # atoms (list): A list of atom lines from a PDB file.
    #
    # Returns:
    # int: The count of hydrogen atoms.
    ####
    hydrogen_count = 0
    for atom in atoms:
        # Extract the atom name which is in columns 13-16 (1-based indexing)
        atom_name = atom[12:16].strip()
        # Check if the atom name starts with "H"
        if atom_name.startswith("H"):
            hydrogen_count += 1
    return hydrogen_count

def verify_protonation(original_pdb, protonated_pdb):
    #### Verifies the protonation by comparing the number of hydrogen atoms 
    # in the original and protonated PDB files.
    #
    # Parameters:
    # original_pdb (str): The path to the original PDB file.
    # protonated_pdb (str): The path to the protonated PDB file.
    #
    # Returns:
    # tuple: A tuple containing the hydrogen counts in the original and protonated PDB files.
    ####
    # Parse the original and protonated PDB files
    original_atoms = parse_pdb(original_pdb)
    protonated_atoms = parse_pdb(protonated_pdb)
    
    # Count the hydrogen atoms in both files
    original_hydrogen_count = count_hydrogens(original_atoms)
    protonated_hydrogen_count = count_hydrogens(protonated_atoms)
    
    return original_hydrogen_count, protonated_hydrogen_count

if __name__ == "__main__":
    # Define the filenames for the original and protonated PDB files
    original_file = "k1_no_H.pdb"
    protonated_file = "###TEMPLATE.pdb"

    # Construct the full paths to the PDB files
    original_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", original_file)
    protonated_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", protonated_file)
    
    # Verify the protonation by comparing hydrogen atom counts
    original_hydrogen_count, protonated_hydrogen_count = verify_protonation(original_path, protonated_path)
    
    # Print the hydrogen atom counts
    print(f"Hydrogen atoms in original file ({original_file}): {original_hydrogen_count}")
    print(f"Hydrogen atoms in protonated file ({protonated_file}): {protonated_hydrogen_count}")
    
    # Check if hydrogen atoms were added in the protonated file
    if protonated_hydrogen_count > original_hydrogen_count:
        print("Verification successful: Hydrogen atoms were added.")
    else:
        print("Verification failed: No hydrogen atoms were added.")
