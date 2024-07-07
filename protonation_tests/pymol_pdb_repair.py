from pymol import cmd

# Load the PDB file
cmd.load("cleaned_k1.pdb", "protein")

# Remove solvent and heteroatoms
cmd.remove("solvent")
cmd.remove("hetatm")

# Add missing hydrogens
cmd.h_add()

# Save the repaired PDB file
cmd.save("repaired_k1.pdb", "protein")
