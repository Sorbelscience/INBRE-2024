import pandas as pd
import pymol
from pymol import cmd
import os
import argparse

# Setup the argument parser
parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, help="directory with .pdb or .cif files", default="/home/storm/INBRE-2024/k1_protonation_states/")
parser.add_argument("--outputDir", type=str, help="name of the directory to save results to (structures and dataset)", default="/home/storm/INBRE-2024/k1_protonation_results/")
parser.add_argument("-c", "--clean", action="store_true", help="if true, clean up files by removing extra chains")
parser.add_argument("-b", "--bsa", action="store_true", help="if true, calculate SASA and BSA")
args = parser.parse_args()

# Data columns for the DataFrame
cols = ['pdb_id', 'chain 1', 'chain 2', 'protein_sasa', 'peptide_sasa', 'combined_sasa', 'bsa']
lst = []

# Create target Directory if don't exist
if not os.path.exists(args.outputDir):
    os.mkdir(args.outputDir)
else:    
    print("Directory " , args.outputDir ,  " already exists")

def pdb_separator():
    """Function to clean PDB files by removing extra chains."""
    filename = "k1.pdb"  # Specifically targeting 'k1.pdb'
    filenamefull = os.path.join(args.inputDir, filename)
    try:
        # Initialize PyMOL
        cmd.reinitialize()
        cmd.load(filenamefull)

        # Get all chains
        chains = cmd.get_chains()
        
        # Clean up unnecessary chains
        for chain in chains:
            if chain != 'A':  # Preserving only chain 'A'
                cmd.remove(f"chain {chain}")
        
        # Save cleaned file
        cleaned_filename = 'clean_' + filename
        cmd.save(os.path.join(args.outputDir, cleaned_filename), "all", -1)
        print("Cleaning complete for ", filename)

    except Exception as e:
        print("Failed at cleaning PDB file", filename, "with error:", str(e))
    return(0)

def calculate_bsa():
    """Function to calculate BSA and SASA."""
    path2 = args.outputDir if args.clean else args.inputDir
    filename = "clean_k1.pdb" if args.clean else "k1.pdb"  # Use cleaned file if cleaned, else original
    filenamefull = os.path.join(path2, filename)

    try:
        # Initialize PyMOL
        cmd.reinitialize()
        cmd.load(filenamefull)

        # Only works if the file actually contains two or more chains
        chains = cmd.get_chains()
        if len(chains) > 1:
            # Create objects for protein, peptide, and protein+peptide pair
            cmd.create("protein", f"{filename} and chain {chains[0]}")
            cmd.create("peptide", f"{filename} and chain {chains[1]}")
            cmd.create("combined", f"{filename} and chain {chains[0]}+{chains[1]}")

            # Add hydrogens and calculate areas
            cmd.h_add("combined")
            protein_area = cmd.get_area("protein")
            peptide_area = cmd.get_area("peptide")
            combined_area = cmd.get_area("combined")

            # Append results
            lst.append([filename, chains[0], chains[1], protein_area, peptide_area, combined_area, (protein_area + peptide_area - combined_area)])
            print("Areas calculated for ", filename)

    except Exception as e:
        print("Failed at calculating BSA for PDB file", filename, "with error:", str(e))
    
    if lst:
        df1 = pd.DataFrame(lst, columns=cols)
        dataset_name = "areas.csv"
        df1.to_csv(os.path.join(args.outputDir, dataset_name))
        print("Data saved to ", dataset_name)

if args.clean:
    pdb_separator()

if args.bsa:
    calculate_bsa()

if not args.bsa and not args.clean:
    print("You need to pass at least one extra argument, -c to clean extra chains and/or -b to calculate SASA and BSA")
