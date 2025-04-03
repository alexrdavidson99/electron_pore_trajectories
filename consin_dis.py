"""
sampling from a cosine distribution
https://www.particleincell.com/2015/cosine-distribution/

@author: lubos brieda
"""
from random import random
import math
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def cosine_dis(xi,r):
    """
    Sample random points from a cosine distribution.
    """
    angle_dis = []
    vel_a = []
    
    # surface properties
    tang1 = [-xi[1]/r, xi[0]/r, 0]
    tang2 = [0, 0, -1]
    norm = [-xi[0]/r, -xi[1]/r, 0]

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


def calculate_theta_cylinder(velocity, impact_position, r):
    
    """
    Calculate the angle θ between the velocity vector of the electron
    and the z-y plane at the moment of impact.
    
    Parameters:
    velocity : ndarray
        The velocity vector of the electron at the moment of impact.
    impact_position : ndarray
        The position of the electron at the moment of impact.
    
    Returns:
    theta : float
        The angle of impact in radians with respect to the z-y plane.
    theta_degrees : float
        The angle of impact in degrees with respect to the z-y plane.
    """
    # Define the radial vector in the z-y plane
    radial_vector = np.array([0, impact_position[1], impact_position[2]])
    
    # Normalize the radial vector
    radial_magnitude = np.linalg.norm(radial_vector)
    radial_normalized = radial_vector / radial_magnitude
    
    # Normalize the velocity vector
    velocity_magnitude = np.linalg.norm(velocity)
    velocity_normalized = velocity / velocity_magnitude
    
    # Compute the dot product between the velocity and the radial vector
    dot_product = np.dot(velocity_normalized, radial_normalized)
    
    # Clip the dot product to avoid out-of-range errors in arccos
    #dot_product = np.clip(dot_product, -1.0, 1.0)
    
    # Calculate the angle with respect to the z-y plane
    theta = np.arccos(dot_product)
    
    # Convert to degrees
    theta_degrees = np.degrees(theta)
    
    return theta, theta_degrees


plot = False

if plot:

    xi = [0.2724217, 0., 0.23779032]
    r = 0.02724217
    v_bin = []
    
    for i in range(0, 10000):
        vel_a, theta = cosine_dis(xi, r)
        

        v_en = 1
        v1 = np.array([vel_a[0][0]*v_en, vel_a[0][1]*v_en, vel_a[0][2]*v_en])
        v_bin.append(v1)
    
    v_bin = np.array(v_bin)
    fig = plt.figure()
    ax2 = fig.add_subplot(212, projection='3d')
    ax2.scatter(v_bin[:, 0],v_bin[:, 1],v_bin[:, 2],c='r',s=0.005);
    ax2.set_xlabel('angle')
    ax2.set_ylabel('normalized count')
    plt.show()


