#!/bin/csh
# Set up environment and variables
set openmm_path = "/path/to/openmm/bin"  # Update this with the actual path to your OpenMM installation
setenv PATH ${openmm_path}:$PATH  # Add OpenMM to the PATH

# Define base directory and input files
set base_dir = "/home/storm/INBRE-2024/charmm-gui-1875490003/openmm"
set init_pdb = "$base_dir/step5_assembly.pdb"
set init_psf = "$base_dir/step5_input.psf"
set init_crd = "$base_dir/step5_input.crd"

# Restraint files
set pos_restraint = "$base_dir/membrane_restraint.str"
set dih_restraint = "$base_dir/membrane_restraint2.str"

# CHARMM Minimization and Setup
set min_str = "$base_dir/step5_input_minimization.str"

# Equilibration and Production
set equi_prefix = "step6"
set prod_inp = "$base_dir/step7_production.inp"

# Run CHARMM Minimization
python -u openmm_run.py -i $min_str -o $base_dir/minimization.out

# Loop through equilibration steps
foreach cnt (1 2 3 4 5 6)
    set equi_inp = "$base_dir/${equi_prefix}.${cnt}_equilibration.inp"
    set equi_out = "$base_dir/${equi_prefix}_${cnt}_equilibration.out"

    # Run equilibration step
    python -u openmm_run.py -i $equi_inp -o $equi_out

    # Check output for errors and stability, adjust as needed
    if (`grep -i "ERROR" $equi_out` != "") then
        echo "Error detected in equilibration step $cnt."
        exit 1
    endif
end

# Run production phase
set prod_out = "$base_dir/production.out"
python -u openmm_run.py -i $prod_inp -o $prod_out

# Check final output
if (`grep -i "ERROR" $prod_out` != "") then
    echo "Error detected in production phase."
    exit 1
endif

echo "Simulation completed successfully."
