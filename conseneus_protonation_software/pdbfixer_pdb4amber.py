import os
import subprocess
import tempfile
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Define the working directory
working_directory = os.path.expanduser("~/INBRE-2024/conseneus_protonation_software")

# Change to the working directory
os.chdir(working_directory)

# Define file paths
input_pdb_file = "k1_propka_protonated.pdb"
output_pdb_file = "k1_protonated_amber_fixed.pdb"

# Use temporary file for the intermediate PDB file
with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, dir=working_directory) as temp_file:
    intermediate_pdb_file = temp_file.name

def add_missing_atoms_and_hydrogens(fixer):
    try:
        fixer.addMissingAtoms()
    except ZeroDivisionError as e:
        print(f"Error adding missing atoms: {e}")

    try:
        fixer.addMissingHydrogens()  # Add missing hydrogens
    except ZeroDivisionError as e:
        print(f"Error adding missing hydrogens: {e}")

try:
    # Step 1: Run PDBFixer to add missing atoms and hydrogens
    fixer = PDBFixer(filename=input_pdb_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    add_missing_atoms_and_hydrogens(fixer)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(intermediate_pdb_file, 'w'))
    print(f"PDBFixer run completed and saved to {intermediate_pdb_file}")

    # Step 2: Run pdb4amber to prepare the PDB file for AMBER
    command = f"pdb4amber -i {intermediate_pdb_file} -o {output_pdb_file} --nohyd"
    subprocess.run(command, shell=True, check=True)
    print(f"pdb4amber run completed and saved to {output_pdb_file}")

finally:
    # Clean up the intermediate file
    os.remove(intermediate_pdb_file)
    print(f"Intermediate file {intermediate_pdb_file} removed")

print(f"Final AMBER-ready PDB file saved as {output_pdb_file}")
