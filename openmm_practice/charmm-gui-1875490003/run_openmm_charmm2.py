from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import os
import time

def main():
    # Define file paths
    psf_path = '/home/storm/INBRE-2024/charmm-gui-1875490003/openmm/step5_input.psf'
    pdb_path = '/home/storm/INBRE-2024/charmm-gui-1875490003/openmm/step5_input.pdb'
    param_files = [
        '/home/storm/INBRE-2024/charmm-gui-1875490003/toppar/par_all36m_prot.prm',
        '/home/storm/INBRE-2024/charmm-gui-1875490003/toppar/top_all36_prot.rtf',
        '/home/storm/INBRE-2024/charmm-gui-1875490003/toppar/par_all36_lipid.prm',
        '/home/storm/INBRE-2024/charmm-gui-1875490003/toppar/top_all36_lipid.rtf',
        '/home/storm/INBRE-2024/charmm-gui-1875490003/toppar/toppar_water_ions.str'
    ]

    # Check if all files exist
    for file_path in [psf_path, pdb_path] + param_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

    # Load CHARMM files
    print("Loading CHARMM PSF file...")
    psf = CharmmPsfFile(psf_path)
    print("Loaded PSF file.")
    
    print("Loading CHARMM PDB file...")
    pdb = PDBFile(pdb_path)
    print("Loaded PDB file.")
    
    print("Loading CHARMM parameters...")
    params = CharmmParameterSet(*param_files)
    print("Loaded CHARMM parameters.")
    
    # Create the system
    print("Creating the system...")
    system = psf.createSystem(params, nonbondedMethod=PME,  # PME for periodic systems
                              nonbondedCutoff=1*nanometer, constraints=HBonds)
    print("System created.")
    
    # Add a barostat for constant pressure
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    
    # Create an integrator
    print("Creating integrator...")
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    print("Integrator created.")
    
    # Create the simulation object
    print("Creating simulation...")
    simulation = Simulation(psf.topology, system, integrator)
    print("Simulation created.")
    
    # Set initial positions
    print("Setting initial positions...")
    simulation.context.setPositions(pdb.positions)
    print("Initial positions set.")
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    print("Energy minimized.")
    
    # Equilibration phases
    equilibration_steps = [
        ('step6.1_equilibration.inp', 50000),
        ('step6.2_equilibration.inp', 50000),
        ('step6.3_equilibration.inp', 50000),
        ('step6.4_equilibration.inp', 50000),
        ('step6.5_equilibration.inp', 50000),
        ('step6.6_equilibration.inp', 50000),
    ]

    for step, nsteps in equilibration_steps:
        print(f"Running equilibration: {step}")
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        simulation.step(nsteps)
        print(f"Equilibration {step} completed.")

    # Set up the reporters
    print("Setting up reporters...")
    simulation.reporters.append(PDBReporter('output.pdb', 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                                                  potentialEnergy=True, temperature=True))
    print("Reporters set up.")
    
    # Run the production simulation
    production_steps = 100000  # Adjust the number of steps as needed
    print("Starting production simulation...")
    simulation.step(production_steps)
    print("Production simulation completed.")

if __name__ == "__main__":
    print("Starting script...")
    main()
