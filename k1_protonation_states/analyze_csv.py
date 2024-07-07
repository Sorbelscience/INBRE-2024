import matplotlib.pyplot as plt
import pandas as pd

# Load the data
data = pd.read_csv('pyPKa_titration.csv', delimiter=';')

# Rename columns if necessary (assuming the columns are 'pH' and 'pKa' based on the previous example)
data.columns = ['pH', 'pKa']

# Set the dark mode for plots
plt.style.use('dark_background')

# Plot 1: pKa vs. pH
plt.figure()
plt.plot(data['pH'], data['pKa'], color='cyan')
plt.xlabel('pH')
plt.ylabel('pKa')
plt.title('pKa vs. pH for pyPKa')
plt.savefig('pyPKa_pKa_vs_pH_dark.png')

# Assuming there are additional data files that need to be processed for other metrics:
# For example, if there are potential energy, temperature, and RMSD data files.

# Load the simulation data from a hypothetical CSV file (this part should be modified based on actual data files)
simulation_data = pd.read_csv('simulation_data.csv')  # Replace with the actual file

# Plot 2: Potential Energy vs. Time
plt.figure()
plt.plot(simulation_data['Time'], simulation_data['PotentialEnergy'], color='cyan')
plt.xlabel('Time (ps)')
plt.ylabel('Potential Energy (kJ/mol)')
plt.title('Potential Energy vs. Time for pyPKa')
plt.savefig('pyPKa_potential_energy_dark.png')

# Plot 3: Temperature vs. Time
plt.figure()
plt.plot(simulation_data['Time'], simulation_data['Temperature'], color='magenta')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Time for pyPKa')
plt.savefig('pyPKa_temperature_dark.png')

# Plot 4: RMSD vs. Time
plt.figure()
plt.plot(simulation_data['Time'], simulation_data['RMSD'], color='yellow')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (nm)')
plt.title('RMSD vs. Time for pyPKa')
plt.savefig('pyPKa_rmsd_dark.png')

plt.show()
