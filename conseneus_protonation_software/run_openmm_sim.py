from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

def run_simulation(top_file, crd_file, output_prefix, num_steps=10000):
    # Load Amber topology and coordinate files
    prmtop = AmberPrmtopFile(top_file)
    inpcrd = AmberInpcrdFile(crd_file)
    
    # Create the system without specifying nonbondedMethod and nonbondedCutoff for non-periodic systems
    system = prmtop.createSystem()
    
    # Create the integrator
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    
    # Set up the simulation
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    
    # Equilibrate the system
    print("Equilibrating...")
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    
    # Set up reporters to record data
    simulation.reporters.append(PDBReporter(f'{output_prefix}_output.pdb', 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, 
                                                  potentialEnergy=True, temperature=True))
    simulation.reporters.append(DCDReporter(f'{output_prefix}_trajectory.dcd', 1000))
    
    # Run the simulation
    print("Running simulation...")
    simulation.step(num_steps)
    
    # Save the final positions
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_final.pdb', 'w'))
    
    print("Simulation completed.")

# Run the simulation
top_file = '0.15_80_10_pH4.6_k1_no_H.top'  # Path to your topology file
crd_file = '0.15_80_10_pH4.6_k1_no_H.crd'  # Path to your coordinate file
output_prefix = 'k1_simulation'  # Prefix for output files
num_steps = 10000  # Number of simulation steps

run_simulation(top_file, crd_file, output_prefix, num_steps)
