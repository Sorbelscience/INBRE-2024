import MDAnalysis as mda
from MDAnalysis.analysis import rms, polymer
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import numpy as np
import re

# Load the modified PDB and DCD files
u = mda.Universe('modified_k1.pdb', 'trajectory.dcd')

# Select all atoms
all_atoms = u.select_atoms('all')

# Parse PROPKA output to get protonation states
def parse_propka_output(filename):
    protonation_states = {}
    with open(filename, 'r') as f:
        for line in f:
            if re.match(r'^\s*\d+\s+[A-Z]{3}\s+[A-Z]\s+\d+\s+\d+\.\d+\s+\d+\.\d+', line):
                parts = line.split()
                residue_name = parts[1]
                chain_id = parts[2]
                residue_number = int(parts[3])
                pka_value = float(parts[4])
                protonation_states[(residue_name, chain_id, residue_number)] = pka_value
    return protonation_states

# RMSD Analysis
def calculate_rmsd():
    R = rms.RMSD(all_atoms, all_atoms, select='all')
    R.run()
    return R.results.rmsd

# Radius of Gyration Analysis
def calculate_radius_of_gyration():
    Rg = []
    for ts in u.trajectory:
        Rg.append(all_atoms.radius_of_gyration())
    return np.array(Rg)

# Hydrogen Bond Analysis
def calculate_hydrogen_bonds():
    donors_sel = 'name N HN'
    acceptors_sel = 'name O'
    hydrogens_sel = 'name H'
    h = HBA(universe=u, donors_sel=donors_sel, acceptors_sel=acceptors_sel, hydrogens_sel=hydrogens_sel)
    h.run()
    return h.count_by_time()

# Plot RMSD
def plot_rmsd(rmsd):
    plt.figure()
    plt.plot(rmsd[:, 1], rmsd[:, 2])
    plt.xlabel('Frame')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD over time')
    plt.savefig('rmsd.png')

# Plot Radius of Gyration
def plot_radius_of_gyration(Rg):
    plt.figure()
    plt.plot(range(len(Rg)), Rg)
    plt.xlabel('Frame')
    plt.ylabel('Radius of Gyration (Å)')
    plt.title('Radius of Gyration over time')
    plt.savefig('radius_of_gyration.png')

# Plot Hydrogen Bonds
def plot_hydrogen_bonds(hbond_counts):
    plt.figure()
    plt.plot(range(len(hbond_counts)), hbond_counts)
    plt.xlabel('Frame')
    plt.ylabel('Number of Hydrogen Bonds')
    plt.title('Hydrogen Bonds over time')
    plt.savefig('hydrogen_bonds.png')

# Highlight Key Residues Based on pKa
def highlight_key_residues(protonation_states, pH=4.6):
    key_residues = {k: v for k, v in protonation_states.items() if abs(v - pH) < 1.0}
    print(f"Key residues near pH {pH}:")
    for res, pka in key_residues.items():
        print(f"Residue: {res}, pKa: {pka}")
    return key_residues

# Main function
def main():
    # Parse PROPKA output
    protonation_states = parse_propka_output('k1.pka')
    
    # Highlight key residues
    key_residues = highlight_key_residues(protonation_states)
    
    # Calculate and plot RMSD
    rmsd = calculate_rmsd()
    plot_rmsd(rmsd)
    
    # Calculate and plot Radius of Gyration
    Rg = calculate_radius_of_gyration()
    plot_radius_of_gyration(Rg)
    
    # Calculate and plot Hydrogen Bonds
    hbond_counts = calculate_hydrogen_bonds()
    plot_hydrogen_bonds(hbond_counts)
    
    print("Analysis complete. Plots saved as 'rmsd.png', 'radius_of_gyration.png', and 'hydrogen_bonds.png'.")

if __name__ == "__main__":
    main()
