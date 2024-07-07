import os
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, stderr, exit
import traceback
import shutil

def main():
    print("Starting script..."); stdout.flush()
    try:
        # Define paths to input files based on your directory listings
        base_dir = os.path.dirname(os.path.abspath(__file__))
        psf_path = os.path.join(base_dir, 'openmm/step5_input.psf')
        pdb_path = os.path.join(base_dir, 'openmm/step5_input.pdb')
        param_files = [
            os.path.join(base_dir, 'toppar/par_all36m_prot.prm'), 
            os.path.join(base_dir, 'toppar/top_all36_prot.rtf'),
            os.path.join(base_dir, 'toppar/par_all36_lipid.prm'), 
            os.path.join(base_dir, 'toppar/top_all36_lipid.rtf'),
            os.path.join(base_dir, 'toppar/toppar_water_ions.str')
        ]

        # Check if all files exist
        print("Checking file paths...")
        print(f"PSF path: {psf_path}")
        print(f"PDB path: {pdb_path}")
        for param in param_files:
            print(f"Param file path: {param}")
        assert os.path.exists(psf_path), f"File not found: {psf_path}"
        assert os.path.exists(pdb_path), f"File not found: {pdb_path}"
        for param in param_files:
            assert os.path.exists(param), f"File not found: {param}"

        print("All files are accessible."); stdout.flush()

        # Load CHARMM files
        print("Loading CHARMM PSF file..."); stdout.flush()
        psf = CharmmPsfFile(psf_path)
        print("Loaded PSF file."); stdout.flush()

        print("Loading CHARMM PDB file..."); stdout.flush()
        pdb = PDBFile(pdb_path)
        print("Loaded PDB file."); stdout.flush()

        print("Loading CHARMM parameters..."); stdout.flush()
        params = CharmmParameterSet(*param_files)
        print("Loaded CHARMM parameters."); stdout.flush()

        # Create an OpenMM system
        print("Creating the system..."); stdout.flush()
        system = psf.createSystem(params, nonbondedMethod=CutoffNonPeriodic,  # CutoffNonPeriodic for non-periodic systems
                                  nonbondedCutoff=1*nanometer, constraints=HBonds)
        print("System created."); stdout.flush()

        print("Creating integrator..."); stdout.flush()
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        print("Integrator created."); stdout.flush()

        print("Creating simulation..."); stdout.flush()
        simulation = Simulation(psf.topology, system, integrator)
        print("Simulation created."); stdout.flush()

        print("Setting initial positions..."); stdout.flush()
        simulation.context.setPositions(pdb.positions)
        print("Initial positions set."); stdout.flush()

        print("Minimizing energy..."); stdout.flush()
        simulation.minimizeEnergy()
        print("Energy minimized."); stdout.flush()

        # Setup reporters to monitor the simulation
        print("Setting up reporters..."); stdout.flush()
        output_dir = os.path.join(base_dir, 'New_Data')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_pdb_path = os.path.join(output_dir, 'output.pdb')
        simulation.reporters.append(PDBReporter(output_pdb_path, 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                                                      potentialEnergy=True, temperature=True))
        print("Reporters set up."); stdout.flush()

        # Run the simulation
        print("Starting simulation..."); stdout.flush()
        simulation.step(10000)
        print("Simulation completed."); stdout.flush()

    except Exception as e:
        print("An error occurred:", file=stderr)
        traceback.print_exc(file=stderr)
        stdout.flush()
        stderr.flush()
        exit(1)

if __name__ == "__main__":
    main()
