import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# Set matplotlib to dark mode
plt.style.use('dark_background')

# Load the modified PDB and DCD files from propka results
u = mda.Universe('k1_propka.pqr', 'propka_trajectory.dcd')

# Select all atoms
all_atoms = u.select_atoms('all')
protein = u.select_atoms('protein')
reference = protein

# Protonation states data (manually provided, updated with propka data)
protonation_states = {
    ('GLU', 'A', 45): 3.84,
    ('NTR', 'A', 45): 8.11,
    ('TYR', 'A', 49): 10.52,
    ('ASP', 'A', 50): 3.26,
    ('LYS', 'A', 51): 10.83,
    ('GLU', 'A', 54): 4.27,
    ('LYS', 'A', 56): None,
    ('ASP', 'A', 57): 3.07,
    ('ASP', 'A', 65): 3.97,
    ('LYS', 'A', 70): 10.58,
    ('ASP', 'A', 81): 4.19,
    ('LYS', 'A', 93): 11.16,
    ('GLU', 'A', 98): 4.08,
    ('LYS', 'A', 100): 11.00,
    ('ASP', 'A', 101): 3.72,
    ('ASP', 'A', 102): 3.59,
    ('ASP', 'A', 106): 2.79,
    ('LYS', 'A', 109): 10.85,
    ('HIS', 'A', 126): 5.88,
    ('HIS', 'A', 127): 6.54,
    ('ASP', 'A', 140): 3.24,
    ('CTR', 'A', 147): 3.34,
    ('NTR', 'B', 235): 7.60,
    ('TYR', 'B', 236): 9.69,
    ('GLU', 'B', 240): 4.01,
    ('HIS', 'B', 241): 7.05,
    ('LYS', 'B', 244): None,
    ('TYR', 'B', 247): 11.46,
    ('ASP', 'B', 253): 3.74,
    ('TYR', 'B', 263): 11.04,
    ('GLU', 'B', 267): 3.94,
    ('LYS', 'B', 268): 11.20,
    ('GLU', 'B', 272): 3.95,
    ('ASP', 'B', 273): 4.40,
    ('GLU', 'B', 274): 4.49,
    ('GLU', 'B', 278): 2.58,
    ('TYR', 'B', 281): 11.69,
    ('TYR', 'B', 282): 10.55,
    ('LYS', 'B', 283): None,
    ('TYR', 'B', 285): None,
    ('LYS', 'B', 296): None,
    ('GLU', 'B', 299): 3.88,
    ('GLU', 'B', 300): 3.22,
    ('ASP', 'B', 303): 3.54,
    ('ASP', 'B', 308): 4.14,
    ('GLU', 'B', 310): 4.22,
    ('ASP', 'B', 313): 2.23,
    ('CTR', 'B', 316): 2.45,
    ('HIS', 'B', 316): 6.57
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
    plt.title('RMSD Over Time for propka')
    plt.legend()
    plt.savefig(f'propka_rmsd_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Radius of Gyration Over Time for propka')
    plt.savefig(f'propka_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Hydrogen Bonds Over Time for propka')
    plt.savefig(f'propka_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('RMSF for propka')
    plt.legend()
    plt.savefig(f'propka_rmsf_dark_{timestamp}.png')
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
    
    print(f"Analysis complete. Plots saved as 'propka_rmsd_dark_{interval_ns}ns_{timestamp}.png', 'propka_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png', 'propka_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png', and 'propka_rmsf_dark_{timestamp}.png'.")

if __name__ == "__main__":
    interval_ns = 1  # Set the desired interval here
    main(interval_ns)
