#!/bin/bash

# Prompt the user for the run name
echo "Enter the run name (it can be anything):"
read RunNumber

# Create main directory for the run
echo "Creating main directory..."
mkdir -p $RunNumber/{ranked_pdb_files,pdb_files}

# Define an array of sequences to process
sequences=("fvqrfssqyv" "viksaikktv" "saikktvshn")

# Define the corresponding target files for each sequence
targets=("no_ligand_k1_rec.trg" "no_ligand_k2_rec.trg" "no_ligand_k1_rec.trg")

# Loop over sequences and process each one
for i in "${!sequences[@]}"; do
  sequence=${sequences[$i]}
  target=${targets[$i]}
  
  # Run ADCP for the current sequence
  echo "Running ADCP on $sequence..."
  adcp -t $target -s $sequence -N 1000 -n 10000000 -o $sequence > $sequence.txt
  
  # Create directory for current sequence and move files
  echo "Sorting files for $sequence..."
  mkdir -p $RunNumber/$sequence
  mv ${sequence}_ranked*.pdb $RunNumber/ranked_pdb_files/
  mv ${sequence}.txt $RunNumber/
  mv ${sequence}*.pdb $RunNumber/pdb_files/
  rm ${sequence}*.out
  
  echo "$sequence processing done."
done

echo "All sequences processed. Done."
