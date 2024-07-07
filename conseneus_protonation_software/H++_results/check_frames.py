import MDAnalysis as mda

# Load the modified PDB and DCD files
u = mda.Universe('k1_H++.pdb', 'H++_trajectory.dcd')

# Print the number of frames
num_frames = u.trajectory.n_frames
print(f"Number of frames in the simulation: {num_frames}")
