from consin_dis import cosine_dis
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import scipy.constants 
from stepping import solve_for_intercept_time, step_position

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.linalg import norm

c = scipy.constants.speed_of_light*1e-6 # in mm/ns
V = 1000   # electrods potential in V
voltage_values = [500, 700, 1000]
d = 700   # pore depth in mm
l_over_d = 46
r = (d/l_over_d)*0.5 

print(f"r = {r} radius in mm")



#m = 511e3
m = 938.272e6

E = V*(c**2)/(d*m)  # electric field acceleration in mm/ns^2
e = 200  # in eV
v = np.sqrt(2*e/m)*c  # velocity in mm/ns
print(1.6e-19 * 2 / 2 * r)
print(f"l/D = {(d/(2*r)):.2f} length to diameter ratio")
orientation = np.array([np.sin(0.13962634), 0., np.cos(0.13962634)])
v0 = v*orientation

orientation = np.array([0., 0., 1])
a0 = E*orientation


# Step particle in field until it reaches boundary or exits
# initial step
x_end = []
y_end = []
z_end = []
angle = []

#energy_values = np.linspace(0.1, 2, num=200)
distance_values = np.linspace(0, d, num=200)
#distance_values = np.arange(0, d + 0.005, 0.005)
#distance_values = distance_values+0.005/2
#energy_values = [0.1, 0.5, 1, 1.5, 2]
count_results = {}
for V in voltage_values:
    count_for_e = []
    E = V*(c**2)/(d*m)  # electric field acceleration in mm/ns^2
    a0 = E*orientation

    for d0 in distance_values:
        number_of_runs = 100000
        count = 0
        for i in range (0, number_of_runs):
            
            hit_pos = np.array([r, 0, d0])
            vel_a, theta = cosine_dis(hit_pos, r)
            emmited_energy = 2 # emmited energy of the ion 
            v = np.sqrt(2*emmited_energy/m)*c
            v1 = np.array([vel_a[0][0] * v, vel_a[0][1] * v, vel_a[0][2] * v])
            angle.append(np.arccos((v1[0]/np.sqrt(v1[0]**2+v1[1]**2))*(-hit_pos[0])/r + (v1[1]/np.sqrt(v1[0]**2+v1[1]**2))*(-hit_pos[1])/r))

            t = solve_for_intercept_time(hit_pos, v1, a0, r)
            end_pos = step_position(hit_pos, v1, a0, t)
            x_end.append(end_pos[0])
            y_end.append(end_pos[1])
            z_end.append(end_pos[2])

            if end_pos[2] >= d:
                count+=1
        count_for_e.append(count/number_of_runs)
    count_results[V] = count_for_e

print(f"Number of particles that exited the pore: {count_for_e}")
end_x = []
end_y = []
end_z = []

colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k', 'w']
for idx, V in enumerate(voltage_values):
 plt.plot((distance_values), count_results[V], color=colors[idx], label=f"V = {V} V")
plt.legend()
plt.xlabel("starting distance from the pore entrance exit (mm)")
plt.ylabel("Probability of exiting the pore")
plt.title("Probability of electrons exiting the pore for different starting distances but keeping energy constant at 1 eV")
plt.yscale('log')

data = {'distance_values': distance_values}

# Add count_results for each voltage to the dictionary
for V in voltage_values:
    data[f'V_{V}'] = count_results[V]

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)
df.to_csv(f'prob_results_{emmited_energy}ev-torch-100k-2ev.csv', index=False)
print(df)
# plt.figure()
# plt.hist(angle, bins=100, alpha=0.75, label=f'angle', density=True,color='r')
# plt.title(f'angle of emmited electrons')
# plt.xlabel('angle')
# plt.ylabel('normalized count')


plt.show()
 

