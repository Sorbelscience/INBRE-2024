import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# Set matplotlib to dark mode
plt.style.use('dark_background')

# Load the modified PDB and DCD files from H++
u = mda.Universe('k1_H++.pdb', 'H++_trajectory.dcd')

# Select all atoms
all_atoms = u.select_atoms('all')
protein = u.select_atoms('protein')
reference = protein

# Protonation states data (manually provided, updated with H++ data)
protonation_states = {
    ('NTGLU', 'A', 1): 8.058,
    ('GLU', 'A', 1): 3.916,
    ('TYR', 'A', 5): 11.545,
    ('ASP', 'A', 6): 3.293,
    ('LYS', 'A', 7): 10.780,
    ('GLU', 'A', 10): 4.574,
    ('LYS', 'A', 12): None,
    ('ASP', 'A', 13): 3.032,
    ('LYS', 'A', 26): 10.497,
    ('ASP', 'A', 37): 4.932,
    ('LYS', 'A', 49): 11.275,
    ('GLU', 'A', 54): 4.096,
    ('ASP', 'A', 57): 3.617,
    ('ASP', 'A', 58): 3.604,
    ('ASP', 'A', 62): 2.543,
    ('LYS', 'A', 65): 10.630,
    ('HIS', 'A', 82): 5.452,
    ('HIS', 'A', 83): 6.103,
    ('ASP', 'A', 96): 3.405,
    ('TYR', 'A', 105): 9.716,
    ('GLU', 'A', 109): 4.116,
    ('HIS', 'A', 110): 7.273,
    ('ASP', 'A', 122): 3.547,
    ('GLU', 'A', 136): 3.944,
    ('LYS', 'A', 137): 11.263,
    ('GLU', 'A', 141): 3.688,
    ('ASP', 'A', 142): 6.921,
    ('GLU', 'A', 143): 4.51,
    ('GLU', 'A', 147): 2.109,
    ('TYR', 'A', 150): None,
    ('TYR', 'A', 151): None,
    ('LYS', 'A', 152): None,
    ('TYR', 'A', 154): None,
    ('LYS', 'A', 165): None,
    ('GLU', 'A', 168): 3.858,
    ('GLU', 'A', 169): 4.846,
    ('ASP', 'A', 172): 3.754,
    ('ASP', 'A', 177): 3.668,
    ('GLU', 'A', 179): 4.256,
    ('ASP', 'A', 182): 1.52,
    ('HIS', 'A', 185): 5.07,
    ('CTHIS', 'A', 185): 1.705
}

# RMSD Analysis
def calculate_rmsd(interval_ns=1):
    # Align the trajectory to the reference structure (first frame)
    align.AlignTraj(u, reference, select='protein and name CA', in_memory=True).run()
    R = rms.RMSD(u, reference, select='protein and name CA')
    R.run()
    return R.results.rmsd[::interval_ns]  # Adjust the interval here

# Radius of Gyration Analysis
def calculate_radius_of_gyration(interval_ns=1):
    Rg = []
    for i, ts in enumerate(u.trajectory):
        if i % interval_ns == 0:
            Rg.append(all_atoms.radius_of_gyration())
    return np.array(Rg)

# Hydrogen Bond Analysis
def calculate_hydrogen_bonds(interval_ns=1):
    donors_sel = 'name N HN'
    acceptors_sel = 'name O'
    hydrogens_sel = 'name H'
    h = HBA(universe=u, donors_sel=donors_sel, acceptors_sel=acceptors_sel, hydrogens_sel=hydrogens_sel)
    h.run()
    return h.count_by_time()[::interval_ns]

# RMSF Analysis
def calculate_rmsf():
    rmsf = RMSF(all_atoms).run()
    return rmsf.rmsf

# Generate a timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Plot RMSD
def plot_rmsd(rmsd, interval_ns):
    print("Plotting RMSD...")
    time = np.arange(0, len(rmsd) * interval_ns, interval_ns)  # Adjust the time axis based on interval
    plt.figure()
    plt.plot(time, rmsd[:, 2], color='cyan')
    avg_rmsd = np.mean(rmsd[:, 2])
    plt.axhline(y=avg_rmsd, color='yellow', linestyle='--', label=f'Avg RMSD: {avg_rmsd:.2f} Å')
    step = max(1, len(time) // 10)  # Ensure the step is at least 1
    for i in range(0, len(time), step):  # Annotate every 10th point or more
        plt.text(time[i], rmsd[i, 2], f'{rmsd[i, 2]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD Over Time for H++')
    plt.legend()
    plt.savefig(f'H++_rmsd_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("RMSD plot saved.")

# Plot Radius of Gyration
def plot_radius_of_gyration(Rg, interval_ns):
    print("Plotting Radius of Gyration...")
    time = np.arange(0, len(Rg) * interval_ns, interval_ns)  # Adjust the time axis based on interval
    plt.figure()
    plt.plot(time, Rg, color='cyan')
    step = max(1, len(time) // 10)  # Ensure the step is at least 1
    for i in range(0, len(time), step):  # Annotate every 10th point or more
        plt.text(time[i], Rg[i], f'{Rg[i]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('Radius of Gyration (Å)')
    plt.title('Radius of Gyration Over Time for H++')
    plt.savefig(f'H++_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("Radius of Gyration plot saved.")

# Plot Hydrogen Bonds
def plot_hydrogen_bonds(hbond_counts, interval_ns):
    print("Plotting Hydrogen Bonds...")
    time = np.arange(0, len(hbond_counts) * interval_ns, interval_ns)  # Adjust the time axis based on interval
    plt.figure()
    plt.plot(time, hbond_counts, color='cyan')
    step = max(1, len(time) // 10)  # Ensure the step is at least 1
    for i in range(0, len(time), step):  # Annotate every 10th point or more
        plt.text(time[i], hbond_counts[i], f'{hbond_counts[i]}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('Number of Hydrogen Bonds')
    plt.title('Hydrogen Bonds Over Time for H++')
    plt.savefig(f'H++_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("Hydrogen Bonds plot saved.")

# Plot RMSF
def plot_rmsf(rmsf):
    print("Plotting RMSF...")
    residues = np.arange(len(rmsf))
    plt.figure()
    plt.plot(residues, rmsf, color='cyan')
    avg_rmsf = np.mean(rmsf)
    plt.axhline(y=avg_rmsf, color='yellow', linestyle='--', label=f'Avg RMSF: {avg_rmsf:.2f} Å')
    step = max(1, len(residues) // 10)  # Ensure the step is at least 1
    for i in range(0, len(residues), step):  # Annotate every 10th residue or more
        plt.text(residues[i], rmsf[i], f'{rmsf[i]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Residue')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF for H++')
    plt.legend()
    plt.savefig(f'H++_rmsf_dark_{timestamp}.png')
    plt.close()
    print("RMSF plot saved.")

# Highlight Key Residues Based on pKa
def highlight_key_residues(protonation_states, pH=4.6):
    key_residues = {k: v for k, v in protonation_states.items() if v and abs(v - pH) < 1.0}
    print(f"Key residues near pH {pH}:")
    for res, pka in key_residues.items():
        print(f"Residue: {res}, pKa: {pka}")
    return key_residues

# Main function
def main(interval_ns=1):
    # Print the number of frames
    num_frames = u.trajectory.n_frames
    print(f"Number of frames in the simulation: {num_frames}")

    # Highlight key residues
    key_residues = highlight_key_residues(protonation_states)
    
    # Calculate and plot RMSD
    rmsd = calculate_rmsd(interval_ns)
    plot_rmsd(rmsd, interval_ns)
    
    # Calculate and plot Radius of Gyration
    Rg = calculate_radius_of_gyration(interval_ns)
    plot_radius_of_gyration(Rg, interval_ns)
    
    # Calculate and plot Hydrogen Bonds
    hbond_counts = calculate_hydrogen_bonds(interval_ns)
    plot_hydrogen_bonds(hbond_counts, interval_ns)
    
    # Calculate and plot RMSF
    rmsf = calculate_rmsf()
    plot_rmsf(rmsf)
    
    print(f"Analysis complete. Plots saved as 'H++_rmsd_dark_{interval_ns}ns_{timestamp}.png', 'H++_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png', 'H++_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png', and 'H++_rmsf_dark_{timestamp}.png'.")

if __name__ == "__main__":
    interval_ns = 1  # Set the desired interval here
    main(interval_ns)

