
from stepping import  step_position, step_energy
import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
from random import random
import math

def cosine_dis():

    angle_dis = []
    vel_a = []
    nbins = 45
    delta_bin = 90 / (nbins)
    # surface properties
    angle_radians = 180
    #angle_radians = 0
    tang1 = np.array([1, 0, 0])  # Tangential direction 1
    tang2 = np.array([0, 0, 1])  # Tangential direction 2
    norm = np.array([0, 1, 0])   # Normal to the y-axis
    # sample random points
    for it in range(0, 1):
        sin_theta = math.sqrt(random())
        cos_theta = math.sqrt(1 - sin_theta * sin_theta)

        # random in plane angle
        psi = random() * 2 * math.pi;

        # three vector components
        a = sin_theta * math.cos(psi)
        b = sin_theta * math.sin(psi)
        c = cos_theta

        angle_dis.append(a)

        # multiply by corresponding directions
        v1 = [a * n for n in tang1]
        v2 = [b * n for n in tang2]
        v3 = [c * n for n in norm]

        # add up to get velocity, vel=v1+v2+v3

        vel = []
        for i in range(0, 3):
            vel.append(v1[i] + v2[i] + v3[i])

        vel_a.append(vel)

        #print(np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2))
        theta = math.acos(cos_theta) * 180 / math.pi;

    return vel_a, theta


def solve_for_intercept_time(x0, v0, acc, target_distance):
    """
    Solve for the intercept time when the particle reaches the target distance in the y direction.
    """
    # Define the polynomial coefficients for the equation of motion in the y direction
    coeffs = [
        0.5 * acc[1],  # t^2 term
        v0[1],  # t term
        x0[1] - target_distance  # constant term
    ]

    # Solve the quadratic equation
    roots = np.roots(coeffs)

    # Filter out complex roots and negative times
    real_roots = roots[np.isreal(roots) & (roots >= 0)]

    if len(real_roots) == 0:
        raise ValueError("No valid intercept time found.")

    # Return the smallest positive real root
    return np.min(real_roots)



angle = []


d_list = [100,200,700]
for j in range(len(d_list)):
    hit_position_x = []
    hit_position_y = []
    hit_time = []
    for i in range(200):
        c = scipy.constants.speed_of_light*1e-6 # in mm/ns
        V = d_list[j]   # electrods potential in V
        d = 1.55  # pore depth in mm

        m = 511e3    # in eV

        E = V*(c**2)/(d*m)  # electric field acceleration in mm/ns^2
        e = 1 # in eV
        v = np.sqrt(2*e/m)*c  # velocity in mm/ns

        orientation = np.array([0, 1,0])
        x0 = np.array([0, 0, 0])
        v0 = v*orientation
        print(v0)
        a0 = E*orientation

        vel_a, theta = cosine_dis()
        new_velocity = np.array([vel_a[0][0] * v, vel_a[0][1] * v, vel_a[0][2] * v])
        # Calculate the energy of the new velocity
        new_energy = 0.5 * m * np.linalg.norm(new_velocity) ** 2 / (c ** 2)
        print(f"Initial energy: {e} eV")
        print(f"New energy: {new_energy} eV")
        print(f"Initial velocity: {v} mm/ns")
        print(f"New velocity: {(new_velocity)} mm/ns")

        print(new_velocity)
        

        ti = solve_for_intercept_time(x0, new_velocity, a0, d)
        print(f"step energy {step_energy(new_velocity, a0, ti, m)}")
        print(f" time to hit {ti}")  
        xi = step_position(x0, new_velocity, a0, ti)
        hit_position_x.append(xi[0])
        hit_position_y.append(xi[2])
        hit_time.append(ti)
    #plt.scatter(hit_position_x, hit_position_y, alpha=0.5, label=f"d = {d_list[j]}")
    bin_edges = np.arange(min(hit_position_x), max(hit_position_x) + 0.05, 0.05)
    #plt.hist(hit_position_x, bins=bin_edges, density=True, alpha=0.6, label=f"d = {d_list[j]}")
    plt.hist(hit_time, bins=30, density=True, alpha=0.6, label=f"V = {d_list[j]} ")
    #plt.title("Impact position")
    #plt.xlabel("position (mm)")
    plt.xlabel("time (ns)")
    plt.ylabel("counts")
plt.legend()
plt.show()