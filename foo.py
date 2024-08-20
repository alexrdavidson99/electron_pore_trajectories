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
d = 46e-2   # pore depth in mm
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
#energy_values = [0.1, 0.5, 1, 1.5, 2]
count_results = {}
for V in voltage_values:
    count_for_e = []
    E = V*(c**2)/(d*m)  # electric field acceleration in mm/ns^2
    a0 = E*orientation

    for d0 in distance_values:
        
        count = 0
        for i in range (0, 10000):
            x1 = np.array([r, 0, d0])
            vel_a, theta = cosine_dis(x1, r)
        
            v = np.sqrt(2*1/m)*c
            v1 = np.array([vel_a[0][0] * v, vel_a[0][1] * v, vel_a[0][2] * v])
            
            
            angle.append(np.arccos((v1[0]/np.sqrt(v1[0]**2+v1[1]**2))*(-x1[0])/r + (v1[1]/np.sqrt(v1[0]**2+v1[1]**2))*(-x1[1])/r))
            t = solve_for_intercept_time(x1, v1, a0, r)
        

            x2 = step_position(x1, v1, a0, t)
            x_end.append(x2[0])
            y_end.append(x2[1])
            z_end.append(x2[2])

            if x2[2] >= d:
                count+=1
        count_for_e.append(count/10000)
    count_results[V] = count_for_e

print(f"Number of particles that exited the pore: {count_for_e}")
end_x = []
end_y = []
end_z = []

colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k', 'w']
for idx, V in enumerate(voltage_values):
    plt.plot(distance_values, count_results[V], color=colors[idx], label=f"V = {V} V")
plt.legend()
plt.xlabel("starting distance from the pore entrance exit (mm)")
plt.ylabel("Probability of exiting the pore")
plt.title("Probability of exiting the pore for different starting distances but keeping energy constant at 1 eV")
plt.yscale('log')

plt.figure()
#plt.plot(energy_values, count_for_e)
plt.xlabel("Energy (eV)")
plt.ylabel("Probability of exiting the pore")
plt.title("Probability of exiting the pore for different voltages")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x_end, z_end, y_end, c='r', marker='o')


# Cylinder parameters
origin = np.array([0, 0, 0])
p0 = np.array([0, 1e-3, 0])
p1 = np.array([0, d, 0])
R = r  # Radius of the cylinder

# Vector in the direction of the axis
v = p1 - p0
mag = norm(v)
v = v / mag

# Two perpendicular vectors to the axis vector
not_v = np.array([1, 0, 0])
if np.allclose(v, not_v):
    not_v = np.array([0, 1, 0])

n1 = np.cross(v, not_v)
n1 /= norm(n1)
n2 = np.cross(v, n1)

# Surface ranges
theta = np.linspace(0, 2 * np.pi, 100)
z = np.linspace(0, mag, 100)
theta, z = np.meshgrid(theta, z)

# Generate coordinates for surface
X = p0[0] + R * np.sin(theta) * n1[0] + R * np.cos(theta) * n2[0] + v[0] * z
Y = p0[1] + R * np.sin(theta) * n1[1] + R * np.cos(theta) * n2[1] + v[1] * z
Z = p0[2] + R * np.sin(theta) * n1[2] + R * np.cos(theta) * n2[2] + v[2] * z

colors = np.where(((Y >= 0) & (Y <= 1e-2)) | ((Y >= 1.8e-1) & (Y <= 2e-1)), 'red', 'aqua')

# Plot the surface
ax.plot_surface(X, Y, Z, facecolors=colors, alpha=0.02)


plt.figure()
x = np.linspace(0, 0.5 * np.pi, 1000)
y = np.cos(x)
plt.plot(x, y)
plt.hist(angle, bins=10, density=True, alpha=0.6, color='g')
plt.show()
 

