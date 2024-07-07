#!/bin/bash

# Define input and output files
input_file="k1.pdb"
output_file="fixed_k1.pdb"
intermediate_file="intermediate_k1.pdb"

# Run PDBFixer to add missing atoms and save the output
pdbfixer $input_file \
  --add-atoms=all \
  --keep-heterogens=all \
  --replace-nonstandard \
  --add-residues \
  --ph=7.0 \
  --positive-ion=Na+ \
  --negative-ion=Cl- \
  --ionic-strength=0.15 \
  --output=$intermediate_file \
  --verbose

# Add TER lines between chains and fix any nonstandard residues manually
python - << EOF
import re

# Read the intermediate file
with open("$intermediate_file", 'r') as infile:
    lines = infile.readlines()

# Write to the final output file
with open("$output_file", 'w') as outfile:
    current_chain = None
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain_id = line[21]
            if current_chain is None:
                current_chain = chain_id
            elif chain_id != current_chain:
                outfile.write('TER\n')
                current_chain = chain_id
        outfile.write(line)
    outfile.write('TER\nEND\n')

# Function to fix nonstandard residue names
def fix_nonstandard_residues(line):
    replacements = {
        'G ': 'DG',
        'A ': 'DA',
        'C ': 'DC',
        'T ': 'DT',
    }
    for old, new in replacements.items():
        line = re.sub(old, new, line)
    return line

# Read the output file again to fix nonstandard residues
with open("$output_file", 'r') as infile:
    lines = infile.readlines()

with open("$output_file", 'w') as outfile:
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            line = fix_nonstandard_residues(line)
        outfile.write(line)
EOF

# Clean up intermediate files
rm $intermediate_file

echo "PDBFixer processing complete. Output written to $output_file."
