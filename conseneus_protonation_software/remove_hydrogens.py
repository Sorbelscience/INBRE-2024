import subprocess
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

input_file = "k1.pdb"
intermediate_file_fixed = "k1_fixed.pdb"
output_file_no_H = "k1_no_H.pdb"

def run_pdbfixer(input_pdb, output_pdb):
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # Add hydrogens assuming pH 7.0 for completeness
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
    print(f"PDBFixer run completed and saved to {output_pdb}")

def run_pdb4amber(input_pdb, output_pdb):
    command = f"pdb4amber -i {input_pdb} -o {output_pdb} --nohyd"
    subprocess.run(command, shell=True, check=True)
    print(f"pdb4amber run completed and saved to {output_pdb}")

if __name__ == "__main__":
    input_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", input_file)
    intermediate_path_fixed = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", intermediate_file_fixed)
    output_path = os.path.join(os.getenv("HOME"), "INBRE-2024", "conseneus_protonation_software", output_file_no_H)
    
    # Run PDBFixer
    run_pdbfixer(input_path, intermediate_path_fixed)
    
    # Run pdb4amber to remove hydrogens
    run_pdb4amber(intermediate_path_fixed, output_path)
