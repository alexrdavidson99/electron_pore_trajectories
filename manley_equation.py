import numpy as np
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use(mplhep.style.LHCb2)
import os

# Constants
K = 0.017        # Arbitrary constant
V = 1       # Fixed initial energy of secondary electron (eV)
alpha = 50   # Length-to-diameter ratio
plt.figure(figsize=(13, 8))
for alpha in [20, 30,40, 50,60]:
    
# Range of total applied voltages (V0)
    V0 = np.linspace(1, 6000, 500)
    v0 = 2000
    # Gain equation
    G = (K * V0**2 / (4 * V * alpha**2)) ** ((4 * V * alpha**2) / V0)
    
    # Plotting
    
    plt.plot(V0, G, label=rf' $\alpha$ = {alpha}')


plt.xlabel('Total Applied Potential Difference $V_0$ (V)')
plt.ylabel('Gain G')

plt.ylim(1, 1e6)
#plt.title('Gain vs Total Applied Voltage')
plt.yscale('log')  # Log scale to better visualize growth
plt.grid(True)
plt.legend(fontsize=26)
plt.tight_layout()

out_fn = r"c:\Users\lexda\PycharmProjects\electron_pore_trajectories\manley_equation_gain_vs_voltage.png"
plt.savefig(out_fn, dpi=300, bbox_inches='tight')
print("Saved figure to:", os.path.abspath(out_fn))
plt.show()
