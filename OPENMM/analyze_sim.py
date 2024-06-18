import numpy as np
import matplotlib.pyplot as plt

# Load the log file data
data = np.loadtxt("md_log.txt", delimiter=',')

# Extract the columns
step = data[:,0]
potential_energy = data[:,1]
temperature = data[:,2]
volume = data[:,3]

# Use the 'dark_background' style to invert colors
plt.style.use('dark_background')

# Plot Potential Energy
plt.plot(step, potential_energy, color='cyan')
plt.xlabel("Step")
plt.ylabel("Potential energy (kJ/mol)")
plt.title("Potential Energy vs. Steps")
plt.show()

# Plot Temperature
plt.plot(step, temperature, color='magenta')
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.title("Temperature vs. Steps")
plt.show()

# Plot Box Volume
plt.plot(step, volume, color='yellow')
plt.xlabel("Step")
plt.ylabel("Volume (nm^3)")
plt.title("Volume vs. Steps")
plt.show()
