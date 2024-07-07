import re

# Define input and output files
input_file = "fixed_k1.pdb"
output_file = "corrected_fixed_k1.pdb"

# Read the input PDB file
with open(input_file, 'r') as infile:
    lines = infile.readlines()

# Function to fix nonstandard residue names and hydrogen atoms
def fix_residues_and_hydrogens(line):
    replacements = {
        'NF': 'N',  # Example: Correct NF to N
        'CY': 'C',  # Example: Correct CY to C
        'NY': 'N',  # Example: Correct NY to N
        'HE2': 'HE2',  # Adjust as necessary
        'HD2': 'HD2',  # Adjust as necessary
        # Add more replacements if needed
    }
    for old, new in replacements.items():
        line = re.sub(rf'\b{old}\b', new, line)
    return line

# Add TER lines between chains and fix residues and hydrogens
with open(output_file, 'w') as outfile:
    current_chain = None
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain_id = line[21]
            if current_chain is None:
                current_chain = chain_id
            elif chain_id != current_chain:
                outfile.write('TER\n')
                current_chain = chain_id
            line = fix_residues_and_hydrogens(line)
        outfile.write(line)
    outfile.write('TER\nEND\n')

print(f"Corrected PDB file written to {output_file}")
