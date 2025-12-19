import numpy as np
import random
import math
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
from scipy.linalg import norm
from gain_sim import run_simulation


gain_spread_mean_list = []
end_positions_z_list = []
range_of_r = [0.003765]
range_of_v = [600]
angles = [0.7] #list(range(2, 40, 5))
number_of_runs = 1

for v in range_of_v:
    for r in range_of_r:
        end_position_threshold = r*2*71.83
        starting_energy = 463    #188.06,13
        starting_angle = angles[0]
        gain_spread = []
        for j in range(number_of_runs):
            count_above_threshold,  end_positions_z = run_simulation(r,v,starting_energy,starting_angle,single_run=False)
            gain_spread.append(count_above_threshold)
            end_positions_z_filtered = [pos for pos in end_positions_z if pos <= end_position_threshold]
            end_positions_z_list.extend(end_positions_z_filtered)

        gain_spread_mean = np.mean(gain_spread)   
        gain_spread_mean_list.append((v, r, gain_spread_mean))
        gain_spread_df = pd.DataFrame(gain_spread, columns=["gain"])
        gain_spread_df.to_csv(f"gain_spread_{v}_{r}.csv", index=False)

        plt.hist(gain_spread, bins=25, alpha=0.75, label=f" samples = {number_of_runs} mean = {gain_spread_mean}",log=True)
        plt.xlabel('Gain (Number of Electrons leaving one pore)')
        plt.ylabel('counts')
        plt.legend()
plt.figure()      

plt.hist(end_positions_z_list,bins=50, alpha=0.75, label=f" samples = {number_of_runs} mean = {gain_spread_mean}",log=True)
print(f"gain spread mean {gain_spread_mean}")

plt.xlim(0,end_position_threshold)
plt.xlabel('z position (mm)')
plt.ylabel('counts')
plt.legend()

print(f"gain spread {gain_spread_mean_list}")
single_run = False
if single_run == False:  
    plt.figure(figsize=(12, 6))

    for v in range_of_v:
        r_values = np.array( [r for v_val, r, gain in gain_spread_mean_list if v_val == v])
        gain_values = [gain for v_val, r, gain in gain_spread_mean_list if v_val == v]
        plt.scatter((0.00165*2*71.83)/(r_values*2), gain_values, label=f'v = {v:.2f}')

    plt.ylabel('Number of Electrons leaving the pore')
    plt.xlabel('L/D')
    plt.legend()
    plt.yscale('log')
plt.show()