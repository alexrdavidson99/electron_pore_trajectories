import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
from scipy.linalg import norm

from energy_dis import accept_reject
from stepping import electron_yield, solve_for_intercept_time, step_position, step_energy
from consin_dis import cosine_dis

from dataclasses import dataclass, field
from typing import List, Optional, Tuple


first_bounce = True
max_electrons = 30000  # Maximum number of electrons to be generated
energy_threshold = 1  # Energy threshold below which electrons are not tracked

@dataclass
class Electron:
    position: List[float]
    impact_energy : float
    angle: float
    bounce_id: int
    initiating_electron: Optional['Electron'] = None
    secondaries: List['Electron'] = field(default_factory=list)
    end_point: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    # To store the end point

    def parent(self):
        return self.initiating_electron

    def gen_secondaries(self, end_point, impact_energy):
        # Generate secondary electrons
        secondaries = []

        poisson_mean = 3 * (impact_energy / 110) * np.exp(1 - (impact_energy / 110))

        pos_0 = np.random.poisson(poisson_mean, 100)
        electrons_yield = np.random.choice(pos_0)
        if e.bounce_id == 0:
            electrons_yield = 3

        for _ in range(electrons_yield):  # Assume each bounce generates two secondary electrons
            secondary = Electron(
                position=end_point,
                impact_energy=self.impact_energy,  # Example impact_energy reduction
                angle=self.angle,  # Example angle change
                bounce_id=self.bounce_id + 1,
                initiating_electron=self
            )
            secondaries.append(secondary)
        return secondaries


def track_electron(electron, r=0.025, m=511e3, V = 1000):

    c = scipy.constants.speed_of_light * 1e-6
    x0 = electron.position  # Initial position
    d = 2 * r * 60
    E = V * (c ** 2) / (d * m)
    field_orientation = np.array([0., 0., 1])
    a0 = E * field_orientation

    if electron.bounce_id == 0:
        starting_electron_energy = electron.impact_energy
        start_orientation = np.array([np.sin(0.13962634), 0., np.cos(0.13962634)])
        v = np.sqrt(2 * starting_electron_energy / m) * c
        v0 = v * start_orientation
        t = solve_for_intercept_time(x0, v0, a0, r)
        hit_position = step_position(x0, v0, a0, t)
        print(f"hit position {hit_position}")
        impact_energy = step_energy(v, a0, t, m)

    else:
        v_initial = np.sqrt((5 / 511e3) * 2 * (c ** 2))
        vel_a, theta = cosine_dis(electron.position , r)
        electron.angle = theta
        v1 = np.array([vel_a[0][0] * v_initial, vel_a[0][1] * v_initial, vel_a[0][2] * v_initial])
        ti = solve_for_intercept_time(electron.position , v1, a0, r)
        print(f" time to hit {ti}")
        hit_position = step_position(electron.position, v1, a0, ti)
        print(f"hit position {hit_position}")
        impact_energy = step_energy(v1, a0, ti, m)


    print(f"impact energy {impact_energy}")

    # Track electron and determine its end point
    # Simplified example: calculate the end point based on start position and angle
    end_point = hit_position
    electron.end_point = end_point
    electron.impact_energy = impact_energy
    return  end_point,impact_energy


def generate_secondary(end_point, impact_energy, electron):
    return electron.gen_secondaries(end_point, impact_energy)


def track_electron_and_generate_secondaries(electron):
    end_point,impact_energy = track_electron(electron)
    if end_point[2] > end_position_threshold :
        return []
    secondaries = generate_secondary(end_point, impact_energy, electron)
    electron.secondaries = secondaries
    return secondaries


def gen_electron():
    # Generate the first electron
    return Electron(position=[0, 0, 0], impact_energy=100, angle=45, bounce_id=0)

electrons = []
if first_bounce:
    first_electron = gen_electron()
    electrons.append(first_electron)

all_electrons = []
total_electron_count = 0
start_positions = []
end_positions = []
end_position_threshold = 0.025*2*60

while electrons:
    e = electrons.pop()

    if e.end_point[2] > end_position_threshold :
        continue
    if e.position[2] > end_position_threshold :
        continue

    new_electrons = track_electron_and_generate_secondaries(e)
    total_electron_count += len(new_electrons)

    if e.bounce_id > 100 :
        break
    if total_electron_count >= max_electrons:
        break  # Stop the simulation if the maximum number of electrons is exceeded

    all_electrons.extend(new_electrons)
    electrons.extend(new_electrons)
    start_positions.append(e.position)
    end_positions.append(e.end_point)

print(f"Total electrons generated: {total_electron_count}")

# Prepare data for plotting impact_energy vs. bounces
bounce_counts = []
energies = []

for electron in all_electrons:
    position = electron.position
    end = electron.end_point
    print(f" start {position}")
    print(f" end {end}")

    bounce_counts.append(electron.bounce_id)
    energies.append(electron.impact_energy)

# Plotting impact_energy vs. bounces


bounce_counts = {}
for electron in all_electrons:
    bounce_id = electron.bounce_id
    print(f"bounce id {bounce_id}")
    if bounce_id in bounce_counts:
        bounce_counts[bounce_id] += 1
    else:
        bounce_counts[bounce_id] = 1

# Prepare data for plotting
bounces = list(bounce_counts.keys())
electron_counts = list(bounce_counts.values())
fig = plt.figure(figsize=(16, 8))  # Increase the size of the figure
#gs = fig.add_gridspec(1, 2, width_ratios=[5, 1])  # Allocate more space to the 3D subplot
#ax = fig.add_subplot(gs[0], projection='3d')
ax = fig.add_subplot(111, projection='3d')

start_positions = np.array(start_positions)
end_positions = np.array(end_positions)

# Plot start positions
ax.scatter(start_positions[:, 0], start_positions[:, 2], start_positions[:, 1], color='blue', label='Start Positions')
# Plot end positions
ax.scatter(end_positions[:, 0], end_positions[:, 2], end_positions[:, 1], color='red', label='End Positions')


ax.set_xlabel('X')
ax.set_ylabel('z')
ax.set_zlabel('y')
ax.set_ylim(0, 0.025*2*60)
ax.set_title('Electron Start and End Positions')
ax.legend()

# Cylinder parameters
origin = np.array([0, 0, 0])
p0 = np.array([0, 1e-3, 0])
p1 = np.array([0, 0.025*2*60, 0])
R = 0.025  # Radius of the cylinder

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
ax.set_box_aspect([1, 5, 1])  # Aspect ratio is 1:1:2

print(f"electron counts {len(end_positions)}")
plt.figure(figsize=(12, 6))
plt.bar(bounces, electron_counts, edgecolor='w')
plt.xlabel('Bounce Count')
plt.ylabel('Number of Electrons')
plt.title('Total Number of Electrons per Bounce')
plt.grid(True)
plt.show()

total_e = []
electron_list = []
SEY= []

electron_counter = 0 # To keep track of unique electron numbers
initial_energy = 0
initial_energy_distribution = accept_reject(1000)






