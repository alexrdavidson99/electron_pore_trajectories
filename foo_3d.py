import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mticker
import numpy as np

# My axis should display 10⁻¹ but you can switch to e-notation 1.00e+01
def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"  # remove int() if you don't use MaxNLocator
    # return f"{10**val:.2e}"      # e-Notation


# Read the CSV files
df_1ev = pd.read_csv('prob_results_1ev.csv')
df_2ev = pd.read_csv('prob_results_2ev.csv')
df_3ev = pd.read_csv('prob_results_3ev.csv')

# Extract distance values and count results
distance_values = 46e-2-df_1ev['distance_values']
count_results_1ev = np.log10(df_1ev.drop(columns=['distance_values']).values)
count_results_2ev = np.log10(df_2ev.drop(columns=['distance_values']).values)
count_results_3ev = np.log10(df_3ev.drop(columns=['distance_values']).values)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot data for 1 eV
for i, voltage in enumerate(df_1ev.columns[1:]):
    ax.plot(distance_values, count_results_1ev[:, i], zs=1, zdir='z', label=f'1 eV, {voltage}')

# Plot data for 2 eV
for i, voltage in enumerate(df_2ev.columns[1:]):
    ax.plot(distance_values, count_results_2ev[:, i], zs=2, zdir='z', label=f'2 eV, {voltage}')

# Plot data for 3 eV
for i, voltage in enumerate(df_3ev.columns[1:]):
    ax.plot(distance_values, count_results_3ev[:, i], zs=3, zdir='z', label=f'3 eV, {voltage}')
    
# Set labels
ax.set_xlabel('Distance Values (mm)')
ax.set_ylabel('Probability of exiting the pore')

ax.set_zlabel('Energy (eV)')
ax.set_title('Ion escape probability for Different Energies and Voltages')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))




# Show legend
ax.legend()

# Show plot
plt.show()